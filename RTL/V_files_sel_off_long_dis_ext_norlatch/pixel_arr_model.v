`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// model for pixel taxel logic + column pulse generator 
//
module pixel_arr_model(
	input rst,
	input [191:0] spk_in_on,
	input [191:0] spk_in_off,	

	input [11:0] s,
	input [15:0] cix_on,
	input [15:0] cix_off,	
	
	output wire [11:0] n_p,
	output wire [15:0] n_cox_on,
	output wire [15:0] n_cox_off,
	
	output wire [15:0] cix_on_del,
	output wire [15:0] cix_off_del	
    );

	wire [191:0] q_latch_on;
	wire [191:0] q_latch_off;	
	
	//wire [9:0] cix_del;
	wire [15:0] cix_on_pulse;
	wire [15:0] cix_off_pulse;
	wire [15:0] cix_pulse;		
	//assign cix_del = cix_on_del | cix_off_del;
	assign cix_on_pulse = cix_on & ~cix_on_del; 
	assign cix_off_pulse = cix_off & ~cix_off_del; 	
	assign cix_pulse = cix_on_pulse| cix_off_pulse;	

	genvar x;
	genvar y;

	generate

		for (y=0; y<12; y=y+1) begin
			for (x=0; x<16; x=x+1) begin
				pixel_model u_on(
					.rst(rst),
					.spk_in(spk_in_on[y*16+x]),

					.acky(s[y]),
					.ackx_pulse(cix_pulse[x]),

					.q_latch(q_latch_on[y*16+x])
				 );	
				 
				 pixel_model u_off(
					.rst(rst),
					.spk_in(spk_in_off[y*16+x]),

					.acky(s[y]),
					.ackx_pulse(cix_pulse[x]),

					.q_latch(q_latch_off[y*16+x])
				 );
				
			end

			
			//assign n_p[x] = ~(q_latch[x*10] | q_latch[x*10+1] | q_latch[x*10+2] | q_latch[x*10+3] | q_latch[x*10+4] | q_latch[x*10+5] | q_latch[x*10+6] | q_latch[x*10+7] | q_latch[x*10+8] | q_latch[x*10+9]);	
			// model default pull-up for undefined inputs
			pd_logic n_p_logic (
				.pd(q_latch_on[y*16] | q_latch_on[y*16+1] | q_latch_on[y*16+2] | q_latch_on[y*16+3] | q_latch_on[y*16+4] | q_latch_on[y*16+5] | q_latch_on[y*16+6] | q_latch_on[y*16+7] | q_latch_on[y*16+8] | q_latch_on[y*16+9] | q_latch_on[y*16+10] | q_latch_on[y*16+11] | q_latch_on[y*16+12] | q_latch_on[y*16+13] | q_latch_on[y*16+14] | q_latch_on[y*16+15] | q_latch_off[y*16] | q_latch_off[y*16+1] | q_latch_off[y*16+2] | q_latch_off[y*16+3] | q_latch_off[y*16+4] | q_latch_off[y*16+5] | q_latch_off[y*16+6] | q_latch_off[y*16+7] | q_latch_off[y*16+8] | q_latch_off[y*16+9] | q_latch_off[y*16+10] | q_latch_off[y*16+11]  | q_latch_off[y*16+12]  | q_latch_off[y*16+13]  | q_latch_off[y*16+14] | q_latch_off[y*16+15]),
				.line(n_p[y])
				);

		end
		
		
		
		for (x=0; x<16; x=x+1) begin					
			//assign n_cox[y] = ~(q_latch[x*10+y] & s[x]);
			// model default pull-up for undefined inputs
			pd_logic n_cox_on_logic(
				.pd((q_latch_on[x] & s[0]) | (q_latch_on[x+16] & s[1]) | (q_latch_on[x+32] & s[2]) | (q_latch_on[x+48] & s[3]) | (q_latch_on[x+64] & s[4])| (q_latch_on[x+80] & s[5])| (q_latch_on[x+96] & s[6])| (q_latch_on[x+112] & s[7]) | (q_latch_on[x+128] & s[8]) | (q_latch_on[x+144] & s[9])  | (q_latch_on[x+160] & s[10])   | (q_latch_on[x+176] & s[11])),
				.line(n_cox_on[x])
				);
				
			pd_logic n_cox_off_logic(
				.pd((q_latch_off[x] & s[0]) | (q_latch_off[x+16] & s[1]) | (q_latch_off[x+32] & s[2]) | (q_latch_off[x+48] & s[3]) | (q_latch_off[x+64] & s[4])| (q_latch_off[x+80] & s[5])| (q_latch_off[x+96] & s[6])| (q_latch_off[x+112] & s[7]) | (q_latch_off[x+128] & s[8]) | (q_latch_off[x+144] & s[9])  | (q_latch_off[x+160] & s[10])   | (q_latch_off[x+176] & s[11])),
				.line(n_cox_off[x])
				);
							
				
			// cix delays for spike step size	
			delay_cell delays_on(
				.cix(cix_on[x]),
				.cix_del(cix_on_del[x])
				);				
				
			// cix delays for spike step size	
			delay_cell delays_off(
				.cix(cix_off[x]),
				.cix_del(cix_off_del[x])
				);					
				
		end
	endgenerate



endmodule
