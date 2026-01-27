`timescale 1ns / 1ps

module address_enc_16(
	input [15:0] ao,   // LSB = bottom (ROW) or right-most (COL) pixel
	output reg [3:0] address
	

    );
	 
	

// LOGIC to generate N-bit address. Assumes only one of the ao bits are asserted/ granted by the arbiter
	always@(*) begin
		if (ao[0] == 1'b1) address = 4'd0;
		else if (ao[1] == 1'b1) address = 4'd1;
		else if (ao[2] == 1'b1) address = 4'd2;		
		else if (ao[3] == 1'b1) address = 4'd3;
		else if (ao[4] == 1'b1) address = 4'd4;
		else if (ao[5] == 1'b1) address = 4'd5;
		else if (ao[6] == 1'b1) address = 4'd6;
		else if (ao[7] == 1'b1) address = 4'd7;
		else if (ao[8] == 1'b1) address = 4'd8;
		else if (ao[9] == 1'b1) address = 4'd9;
		else if (ao[10] == 1'b1) address = 4'd10;
		else if (ao[11] == 1'b1) address = 4'd11;
		else if (ao[12] == 1'b1) address = 4'd12;
		else if (ao[13] == 1'b1) address = 4'd13;
		else if (ao[14] == 1'b1) address = 4'd14;
		else if (ao[15] == 1'b1) address = 4'd15;   
	end

endmodule
