"""
Defines Convolutional / MLP SNN model (SVM-only hidden spike output)
Gesture recognition using 8GHz Radar
Author: Ali Safa - IMEC- KU Leuven, Federico Corradi - IMEC-NL
Modified by: Mark Alea - KU Leuven
Further modified: GPU-safe tensors + torch-only surrogate grad backward.
04/02/2023
"""

import math
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

torch.autograd.set_detect_anomaly(True)

gamma = 0.5  # gradient scale
lens = 0.5
decay_neu = 1.0  # 1 for IF, <1 for LIF


# ---- Surrogate firing function (torch-only backward: GPU-safe) ----
class ActFun(torch.autograd.Function):
    @staticmethod
    def forward(ctx, input):
        ctx.save_for_backward(input)
        return input.gt(0).to(input.dtype)

    @staticmethod
    def backward(ctx, grad_output):
        (input,) = ctx.saved_tensors
        grad_input = grad_output

        # Gaussian surrogate gradient: exp(-x^2/(2*lens^2)) / (sqrt(2*pi)*lens)
        # Torch-only => stays on GPU if input is on GPU.
        denom = (math.sqrt(2.0 * math.pi) * lens)
        temp = torch.exp(-(input * input) / (2.0 * lens * lens)) / denom

        return gamma * grad_input * temp


act_fun = ActFun.apply


def mem_update_adp(ops, x, mem, thr):
    mem = decay_neu * mem + ops(x)
    inputs_ = mem - thr
    spike = act_fun(inputs_)
    mem = mem * (1 - spike)  # reset
    negative_ = (mem < -thr).to(mem.dtype)
    mem = (mem * (1 - negative_)) - (thr * negative_)
    return mem, spike


def mem_update_NU_adp(inputs, mem, thr):
    mem = decay_neu * mem + inputs
    inputs_ = mem - thr
    spike = act_fun(inputs_)
    mem = mem * (1 - spike)  # reset
    negative_ = (mem < -thr).to(mem.dtype)
    mem = (mem * (1 - negative_)) - (thr * negative_)
    return mem, spike


def moving_filt(window_size, stride, hidden_spike_t):
    # hidden_spike_t: [B, H, T]
    windows = hidden_spike_t.unfold(2, window_size, stride)  # [B,H,Tw,window]
    means = windows.mean(dim=3)
    return means


class mini_eMLP(nn.Module):
    def __init__(self, input_ch, input_x, input_y, hidden_size_1, output_size, window_size):
        super(mini_eMLP, self).__init__()

        self.hidden_size_1 = hidden_size_1
        self.output_size = output_size
        self.window_size = window_size

        input_size = input_ch * input_x * input_y
        self.i2h = nn.Linear(input_size, hidden_size_1)
        self.h2o = nn.Linear(hidden_size_1, output_size)
        self.h2h = nn.Linear(hidden_size_1, hidden_size_1)
        self.thr_h_1 = nn.Parameter(torch.Tensor(hidden_size_1))
        self.thr_o = nn.Parameter(torch.Tensor(output_size))

        nn.init.orthogonal_(self.h2h.weight)
        nn.init.xavier_uniform_(self.i2h.weight)
        nn.init.xavier_uniform_(self.h2o.weight)
        nn.init.constant_(self.i2h.bias, 0)
        nn.init.constant_(self.h2h.bias, 0)
        nn.init.constant_(self.h2o.bias, 0)
        nn.init.constant_(self.thr_h_1, 1.0)
        nn.init.constant_(self.thr_o, 1.0)

    def forward(self, input, labels, hidden_mem, hidden_spike, output_mem, output_spike):
        # input: [B, T, input_size]
        device = input.device
        dtype = input.dtype

        batch_size, T_len, _ = input.shape

        nbr_events_avrg = torch.zeros(T_len, device=device, dtype=dtype)
        output_spike_sum = torch.zeros(batch_size, self.output_size, device=device, dtype=dtype)
        hidden_spike_t = torch.zeros(batch_size, self.hidden_size_1, T_len, device=device, dtype=dtype)

        for this_t in range(T_len):
            input_x = input[:, this_t, :]
            h_input = self.i2h(input_x) + self.h2h(hidden_spike)
            hidden_mem, hidden_spike = mem_update_NU_adp(h_input, hidden_mem, self.thr_h_1)
            output_mem, output_spike = mem_update_adp(self.h2o, hidden_spike, output_mem, self.thr_o)

            hidden_spike_t[:, :, this_t] = hidden_spike
            nbr_events_avrg[this_t] = (output_spike.sum() + hidden_spike.sum()) / float(batch_size)
            output_spike_sum = output_spike_sum + output_spike

        hidden_spike_filt = moving_filt(window_size=self.window_size, stride=1, hidden_spike_t=hidden_spike_t)
        return hidden_spike_filt

