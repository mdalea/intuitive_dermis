"""
Defines 1-Hidden Layer SNN model
Texture and Flutter Recognition
Author: Ali Safa - IMEC- KU Leuven, Federico Corradi - IMEC-NL
Modified by: Mark Alea - KU Leuven
Further modified to be CUDA-safe + faster (no NumPy in backward, no CPU sync in forward)
04/02/2023
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

torch.autograd.set_detect_anomaly(True)

gamma = 0.5   # gradient scale
lens = 0.5
decay_neu = 1 # 1 for IF, <1 for LIF

# define approximate firing function
class ActFun(torch.autograd.Function):
    @staticmethod
    def forward(ctx, input):  # input = membrane potential - threshold
        ctx.save_for_backward(input)
        return input.gt(0).to(input.dtype)  # is firing (float tensor)

    @staticmethod
    def backward(ctx, grad_output):  # approximate the gradients
        (input,) = ctx.saved_tensors
        grad_input = grad_output.clone()

        # CUDA-safe Gaussian surrogate gradient (torch ops, not numpy)
        # temp = exp(-(input^2)/(2*lens^2)) / (sqrt(2*pi)*lens)
        lens_t = torch.as_tensor(lens, device=input.device, dtype=input.dtype)
        denom = torch.sqrt(torch.as_tensor(2.0 * np.pi, device=input.device, dtype=input.dtype)) * lens_t
        temp = torch.exp(-(input * input) / (2.0 * lens_t * lens_t)) / denom

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

class mini_eMLP(nn.Module):
    def __init__(self, input_ch, input_x, input_y, hidden_size_1, output_size, criterion, batch_size):
        super(mini_eMLP, self).__init__()

        self.criterion = criterion
        self.batch_size = batch_size
        self.hidden_size_1 = hidden_size_1
        self.output_size = output_size

        input_size = input_ch * input_x * input_y
        self.i2h = nn.Linear(input_size, hidden_size_1)
        self.h2o = nn.Linear(hidden_size_1, output_size)
        self.h2h = nn.Linear(hidden_size_1, hidden_size_1)

        # learn thresholds
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
        """
        input:  [B, T, input_size]
        labels: [B]
        hidden_mem/spike: [hidden]
        output_mem/spike: [out]
        """

        device = input.device
        dtype = input.dtype

        batch_size, T_len, _ = input.shape

        # allocate on same device (important for CUDA speed)
        nbr_events_avrg = torch.zeros(T_len, device=device, dtype=torch.float32)
        output_spike_sum = torch.zeros(batch_size, self.output_size, device=device, dtype=torch.float32)

        # Ensure state tensors are on the right device
        hidden_mem = hidden_mem.to(device=device, dtype=torch.float32)
        hidden_spike = hidden_spike.to(device=device, dtype=torch.float32)
        output_mem = output_mem.to(device=device, dtype=torch.float32)
        output_spike = output_spike.to(device=device, dtype=torch.float32)

        # BPTT over time
        for this_t in range(T_len):
            input_x = input[:, this_t, :].to(torch.float32)

            h_input = self.i2h(input_x) + self.h2h(hidden_spike)
            hidden_mem, hidden_spike = mem_update_NU_adp(h_input, hidden_mem, self.thr_h_1)
            output_mem, output_spike = mem_update_adp(self.h2o, hidden_spike, output_mem, self.thr_o)

            nbr_events_avrg[this_t] = (output_spike.sum() + hidden_spike.sum()) / batch_size
            output_spike_sum = output_spike_sum + output_spike

        # classification
        output_sumspike = F.log_softmax(output_spike_sum, dim=1)
        pred_ = output_sumspike.argmax(dim=1)  # [B] on-device tensor

        # loss on the same device (no .to(device) inside)
        loss_h = self.criterion(output_sumspike, labels.long())

        # Return predictions as tensor to avoid CPU sync; caller can compute accuracy
        # To match your old "predictions_t = torch.tensor(predictions)" shape, we return [1, B]
        predictions_t = pred_.unsqueeze(0)

        return predictions_t, loss_h, nbr_events_avrg

