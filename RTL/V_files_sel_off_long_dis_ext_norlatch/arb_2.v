//Verilog HDL for "spike_enc", "arb_2" "functional"


module arb_2 (

input l1i,
input l2i,
output wire n_l1o,
output wire n_l2o,
output wire ro,
input n_ri


);

wire a1p;
wire n_a1p,n_a2p;
wire _l1o,_l2o;

// safer to use SR latch instead to take into account when both S & R == 0
// nand sr latch at the input
//assign n_a2p = ~(n_a1p & l2i);
//assign n_a1p = ~(n_a2p & l1i);

//sr_latch_v u0 (.set( l1i ),
sr_latch_nor u0 (.set( l1i ),
			.reset( l2i ),
			.q( a1p ) );


assign n_a1p = ~a1p;
assign n_a2p = a1p;


// l1 or l2 ack'd only when arbiter branch is requesting and ack'd while either l1 or l2 issued request 
assign _l1o = ~(n_ri | n_a1p);
assign _l2o = ~(n_ri | n_a2p);

assign n_l1o = ~(ro & _l1o);
assign n_l2o = ~(ro & _l2o);

//assign ro = n_ri & (l1i | l2i);

// safer to use SR latch instead of asymmetric logic as in Boahen. Allows for synthesis tool
// only send REQ upstream if not yet ack'd and at least one request from l1i or l2i
// only reset REQ upstream if both l1i & l2i goes LOW


//sr_latch_v u1 (.set( (l1i | l2i) ),
sr_latch_nor u1 (.set( (l1i | l2i) ),
			.reset( ~(l1i | l2i) ),
			.q( ro) );



endmodule
