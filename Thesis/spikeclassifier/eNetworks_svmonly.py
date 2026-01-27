"""
Defines Convolutional SNN model
Gesture recognition using 8GHz Radar
Author: Ali Safa - IMEC- KU Leuven, Federico Corradi - IMEC-NL
Modified by: Mark Alea - KU Leuven
04/02/2023
"""
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable



device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
torch.autograd.set_detect_anomaly(True)
gamma = .5  # gradient scale
lens = 0.5
decay_neu = 1 #0.5 #1 for IF, <1 for LIF


# define approximate firing function
class ActFun(torch.autograd.Function):
    @staticmethod
    def forward(ctx, input):  # input = membrane potential- threshold
        ctx.save_for_backward(input)
        return input.gt(0).float()  # is firing
    
    @staticmethod
    def backward(ctx, grad_output):  # approximate the gradients
        input, = ctx.saved_tensors
        grad_input = grad_output.clone()
        temp = np.exp((-(input) ** 2 / (2 * lens * lens))) / (np.sqrt(2 * np.pi) * lens)
        return gamma * grad_input * temp.float()

act_fun = ActFun.apply

def mem_update_adp(ops, x, mem, thr ):
    mem = decay_neu*mem + ops(x)
    inputs_ = mem - thr
    spike = act_fun(inputs_)#(inputs_ > 0).float()
    mem = mem*(1-spike) # reset
    negative_ = (mem < -thr).float()
    mem = (mem*(1-negative_))-(thr*negative_)
    return mem, spike

def mem_update_NU_adp(inputs, mem, thr):
    mem = decay_neu*mem + inputs
    inputs_ = mem - thr
    spike = act_fun(inputs_)#(inputs_ > 0).float()
    mem = mem*(1-spike) # reset
    negative_ = (mem < -thr).float()
    mem = (mem*(1-negative_))-(thr*negative_)
    return mem, spike

def moving_filt(window_size, stride, hidden_spike_t):

    windows = hidden_spike_t.unfold(2, window_size, stride)
    means = windows.mean(dim=3)

    return means


       
class mini_eMLP(nn.Module):
    def __init__(self, input_ch, input_x, input_y, hidden_size_1, output_size, window_size):
        super(mini_eMLP, self).__init__()
        
        #self.batch_size = batch_size
        self.hidden_size_1 = hidden_size_1
        self.output_size = output_size
        self.window_size = window_size

        input_size = input_ch * input_x * input_y
        self.i2h = nn.Linear(input_size, hidden_size_1)
        self.h2o = nn.Linear(hidden_size_1, output_size)
        self.h2h = nn.Linear(hidden_size_1, hidden_size_1)
        self.thr_h_1 = nn.Parameter(torch.Tensor(hidden_size_1))#, requires_grad=True) #learn threshold
        self.thr_o = nn.Parameter(torch.Tensor(output_size))

        nn.init.orthogonal_(self.h2h.weight) #initialize weights following Xavier's rule (maintains equal activity variance along the network)
        nn.init.xavier_uniform_(self.i2h.weight)
        nn.init.xavier_uniform_(self.h2o.weight)
        nn.init.constant_(self.i2h.bias, 0)
        nn.init.constant_(self.h2h.bias, 0)
        nn.init.constant_(self.h2o.bias, 0)
        nn.init.constant_(self.thr_h_1, 1.0)
        nn.init.constant_(self.thr_o, 1.0)


    def forward(self, input, labels, hidden_mem, hidden_spike, output_mem, output_spike):

        # Feed in the whole sequence
        batch_size, T_len, _ = input.shape
        nbr_events_avrg = torch.zeros(T_len, requires_grad=False)

        output_spike_sum = torch.zeros(batch_size, self.output_size)
        hidden_spike_t = torch.zeros(batch_size, self.hidden_size_1, T_len, requires_grad=False)

        for this_t in range(T_len): #BPTT
            #sample randomly the input - 
            input_x = input[:, this_t, :] 

            h_input = self.i2h(input_x.float()) + self.h2h(hidden_spike)            
            hidden_mem, hidden_spike = mem_update_NU_adp(h_input,hidden_mem, self.thr_h_1)
            output_mem, output_spike = mem_update_adp(self.h2o, hidden_spike, output_mem, self.thr_o)
            hidden_spike_t[:,:,this_t] = hidden_spike
            nbr_events_avrg[this_t] = (torch.sum(output_spike) + torch.sum(hidden_spike))/batch_size
            output_spike_sum[:,:] = output_spike_sum[:,:] + output_spike

        # Perform moving window averaging on hidden layer spikes
        hidden_spike_filt = moving_filt(window_size=self.window_size, stride=1, hidden_spike_t=hidden_spike_t)

        return hidden_spike_filt


    def num_flat_features(self, x):
        size = x.size()[1:]  # all dimensions except the batch dimension
        num_features = 1
        for s in size:
            num_features *= s
        return num_features





        
    
    
