`timescale 1ns / 1ps

module aer_row_off_long_dis_ext_norlatch(
	
  	//output n_arb_y_ri,  // FOR SYNTHESIS CONSTRAINT

	input RST,
	input ACK,
	
	input GREEDY,
	
	output wire [3:0] ADDRY,
	
	input [11:0] N_P,
	
	output wire [11:0] S,

	//output arb_y_ro
	input AER_DIS,

	input [10:0] S_EXT
    );

wire [11:0] n_ri_y; // row requests / arbiter input
wire [11:0] ro_y; // row acks / arbiter output
wire [11:0] ao_y; // y-address encoder input

wire arb_y_ro; // y-arbiter request
wire n_arb_y_ri; // y-arbiter ack

wire [11:0] s_int;

assign S[10:0] = s_int[10:0] | S_EXT[10:0];
assign S[11] = s_int[11];

rowcol_interface_12 row_interface ( 
			.rst( RST ),
			.greedy( GREEDY ),
			.n_p( N_P ),
			.n_ai ( ~ACK ),  // neg-asserted ai input
			.n_ri ( n_ri_y ),

			.arbtop_n_ri(n_arb_y_ri),	
			
			.ro ( ro_y ),
			.s ( s_int ),
			.ao ( ao_y ),
			.aer_dis ( AER_DIS )
			);
			
address_enc_12 enc_y(
	.ao ( ao_y ),   
	.address ( ADDRY )


    );		

assign n_arb_y_ri = ~arb_y_ro;
//assign n_arb_y_ri = greedy ? 1'b1 : ~arb_y_ro;
arb_12 arb_y ( 

	.lni ( ro_y ),		
	.n_lno ( n_ri_y ),  
	.ro ( arb_y_ro ),
	.n_ri ( n_arb_y_ri )

);	 
			

	 
endmodule
