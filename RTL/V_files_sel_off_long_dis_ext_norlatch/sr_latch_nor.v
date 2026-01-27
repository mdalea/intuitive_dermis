//Verilog HDL for "spike_enc", "sr_latch_v" "functional"


module sr_latch_nor (set,reset,q);
	
input set,reset;
output q;

wire q, qb;

//nor(q, reset, qb);
//nor(qb, set, q);

// as synthesized
assign q = ~((~q  & ~set) | reset);

endmodule
