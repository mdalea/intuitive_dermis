//Verilog HDL for "spike_enc", "sr_latch_v" "functional"


module sr_latch_v (set,reset,q);
	
input set,reset;
output q;

reg q;

always@(*) begin
	if (set==1'b1)
		//q <= #10 1'b1;     //add delay for debug
		q <= 1'b1;
	else if (reset==1'b1)
		//q <= #10 1'b0;
		q <= 1'b0;
	else
		q <= q;
end

endmodule
