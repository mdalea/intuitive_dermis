`timescale 1ns / 1ps

module aer_col_del(
	
  	//output n_arb_x_ri,  // FOR SYNTHESIS CONSTRAINT

	input RST,
	input ACK,
	
	//input GREEDY,
	output wire REQON,
	output wire REQON_B,	
	output wire REQOFF,
	output wire REQOFF_B,	
	
	//output wire [3:0] ADDRX,
	output wire [2:0] ADDRX,
	output wire [2:0] ADDRX_B,

	//input [15:0] N_COX_ON,
	//input [15:0] N_COX_OFF,
	input [7:0] N_COX_ON,
	input [7:0] N_COX_OFF,	

	//output wire [15:0] CIX_ON,
	//output wire [15:0] CIX_OFF,	
	output wire [7:0] CIX_ON,
	output wire [7:0] CIX_OFF,		
	// for differential per-taxel logic	
	output wire [7:0] CIX_ON_B,
	output wire [7:0] CIX_OFF_B,			
	
	input [7:0] CIX_ON_DEL, // delay taxel reset + delay REQ 
	input [7:0] CIX_OFF_DEL, // delay taxel reset + delay REQn
 	
	input AER_DIS
    );


wire [7:0] n_ri_x; // col requests / arbiter input
wire [7:0] ro_x; // col acks / arbiter output
wire [7:0] ao_x_on; // x-address ON encoder input
wire [7:0] ao_x_off; // x-address OFF encoder input
wire [7:0] ao_x; // OR'd ao's

wire arb_x_ro; // x-arbiter request
wire n_arb_x_ri; // x-arbiter ack

//wire [9:0] cix;
//wire [9:0] n_cox;
// DERMIS -delay is retained to slow down AER
wire cix_on_del_or;
wire cix_off_del_or;

assign REQON_B = ~REQON;
assign REQOFF_B = ~REQOFF;
assign ADDRX_B = ~ADDRX;

assign CIX_ON_B = ~CIX_ON;
assign CIX_OFF_B = ~CIX_OFF;

assign cix_on_del_or = CIX_ON_DEL[0] | CIX_ON_DEL[1] | CIX_ON_DEL[2] | CIX_ON_DEL[3] | CIX_ON_DEL[4] | CIX_ON_DEL[5] | CIX_ON_DEL[6] | CIX_ON_DEL[7] ;
assign cix_off_del_or = CIX_OFF_DEL[0] | CIX_OFF_DEL[1] | CIX_OFF_DEL[2] | CIX_OFF_DEL[3] | CIX_OFF_DEL[4] | CIX_OFF_DEL[5] | CIX_OFF_DEL[6] | CIX_OFF_DEL[7];
//assign cix_on_del_or = CIX_ON_DEL[0] | CIX_ON_DEL[1] ;
//assign cix_off_del_or = CIX_OFF_DEL[0] | CIX_OFF_DEL[1] ;


// REQUEST LOGIC -> goes to receiver chip
//sr_latch_v u0 (
sr_latch_nor u0 (
			.set( cix_on_del_or ),   // assert when row/col line is pulled DOWN
			.reset( ~cix_on_del_or ),  // deasert when neg-asserted address encoder ack is received and row/col line is pulled up
			.q( REQON ) );         // request goes to arbiter

//sr_latch_v u1 (
sr_latch_nor u1 (
			.set( cix_off_del_or ),   // assert when row/col line is pulled DOWN
			.reset( ~cix_off_del_or ),  // deasert when neg-asserted address encoder ack is received and row/col line is pulled up
			.q( REQOFF ) );         // request goes to arbiter			

wire [7:0] ro_x_on;
wire [7:0] ro_x_off;
assign ro_x = ro_x_on | ro_x_off;  // combine request of ON and OFF taxels
			
rowcol_interface_8 col_interface_on ( 
			.rst( RST ),
			//.greedy( GREEDY ),			
			.n_p( N_COX_ON ),
			.n_ai ( ~ACK ),
			.n_ri ( n_ri_x ),

			.arbtop_n_ri(n_arb_x_ri),	
			
			.ro ( ro_x_on ),
			.s ( CIX_ON ),
			.ao ( ao_x_on ),
			.aer_dis (AER_DIS)
			);			

rowcol_interface_8 col_interface_off ( 
			.rst( RST ),
			//.greedy( GREEDY ),			
			.n_p( N_COX_OFF ),
			.n_ai ( ~ACK ),
			.n_ri ( n_ri_x ),

			.arbtop_n_ri(n_arb_x_ri),	
			
			.ro ( ro_x_off ),
			.s ( CIX_OFF ),
			.ao ( ao_x_off ),
			.aer_dis (AER_DIS)
			);

assign ao_x = ao_x_on | ao_x_off;
//address_enc_16 enc_x(
address_enc_8 enc_x(
	.ao ( ao_x ),   
	.address ( ADDRX )

    );
	 
assign n_arb_x_ri = ~arb_x_ro;
//assign n_arb_x_ri = greedy ? 1'b1 : ~arb_x_ro;
//arb_16 arb_x ( 
arb_8 arb_x ( 

	.lni ( ro_x ),		
	.n_lno ( n_ri_x ),  
	.ro ( arb_x_ro ),
	.n_ri ( n_arb_x_ri )

);	 	 
	 
endmodule
