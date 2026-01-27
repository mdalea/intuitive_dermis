module col_or(

	input [15:0] ackx_pulse_ext, // use only 9:0

	input [15:0] ackx_pulse_int,

	output wire [15:0] ackx_pulse

);


assign ackx_pulse[15:0] = ackx_pulse_int[15:0] |  ackx_pulse_ext[15:0];
//assign  ackx_pulse[15:10] =  ackx_pulse_int[15:10] | ackx_pulse_ext[15:0];

endmodule

