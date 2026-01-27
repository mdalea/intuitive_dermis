`timescale 1ns / 1ps

module pd_logic(
	input pd,
	output reg line
    );

	initial begin
		line = 1'b1;
		#10000000;
	end

//pullup u0(line);
//
//assign line = pd ? 1'b0 : 1'bz;

	always@(*) begin
		case(pd)
			1'b1 : line = #1 1'b0;  // make the delay weird ;) so it's easy to trace
			default : line = #1 1'b1;
		endcase
	end

endmodule
