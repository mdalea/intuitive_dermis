//Verilog HDL for "spike_enc", "sr_latch_v" "functional"


module sr_latch_nand (set,reset,q);
	
input set,reset;
output q;

wire q, qb;

nand(q, set, qb);
nand(qb, reset, q);

endmodule
