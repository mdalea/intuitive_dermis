//Verilog HDL for "spike_enc", "rowcol_interface" "functional"



module rowcol_interface (

input rst,
input greedy,
input n_p,
input n_ai,
input n_ri,

input arbtop_n_ri,

output wire ro,
output wire s,
output wire ao,

input aer_dis

);
wire arbtop_n_ri_sel;
assign arbtop_n_ri_sel = greedy ? 1'b1 : arbtop_n_ri;  // if 1, assert new taxel request even without waiting for previous

// REQUEST LOGIC -> goes to arbiter
//sr_latch_v u0 (.set( ~aer_dis & ~n_p & n_ai & arbtop_n_ri_sel),   // assert when row/col line is pulled DOWN & arbiter is finished with previous spike & prevent fast arbiter actions when ACK==1
sr_latch_nor u0 (.set( ~aer_dis & ~n_p & n_ai & arbtop_n_ri_sel),
			.reset( rst | (n_p & ~n_ai) ),  // put rst to propagate to the arbiter
			.q( ro ) );         // request goes to arbiter


// ACK LOGIC -> goes to taxel logic
//sr_latch_v u1 (.set( n_ai & ~n_ri & ~n_p ),  // assert when neg-asserted ACK'd is received from arbiter and address encoder is done with previous
sr_latch_nor u1 (.set( n_ai & ~n_ri & ~n_p ),  // assert when neg-asserted ACK'd is received from arbiter and address encoder is done with previous
// added n_p input to differentiate cix_on and cix_off generation from s
			.reset(  n_ri ),			// deassert when arbiter ACK is deasserted 
			.q( s ) );

assign ao = s;

endmodule
