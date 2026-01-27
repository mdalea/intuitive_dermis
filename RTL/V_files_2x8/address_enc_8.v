`timescale 1ns / 1ps

module address_enc_8(
	input [7:0] ao,   // LSB = bottom (ROW) or right-most (COL) pixel
	output reg [2:0] address

    );



// No latch in Yari's pragmatic cells, change style	 
	
/*
// LOGIC to generate N-bit address. Assumes only one of the ao bits are asserted/ granted by the arbiter
	always@(*) begin
		if (ao[0] == 1'b1) address = 2'd0;
		else if (ao[1] == 1'b1) address = 2'd1;
		else if (ao[2] == 1'b1) address = 2'd2;		
		//else if (ao[3] == 1'b1) address = 4'd3;
		//else if (ao[4] == 1'b1) address = 4'd4;
		//else if (ao[5] == 1'b1) address = 4'd5;
		//else if (ao[6] == 1'b1) address = 4'd6;
		//else if (ao[7] == 1'b1) address = 4'd7;
		//else if (ao[8] == 1'b1) address = 4'd8;
		//else if (ao[9] == 1'b1) address = 4'd9;
		//else if (ao[10] == 1'b1) address = 4'd10;
		//else if (ao[11] == 1'b1) address = 4'd11;
		////else address = 4'd15;   // handle erroneous address encoding
	end
*/

/*
// LOGIC to generate N-bit address. Assumes only one of the ao bits are asserted/ granted by the arbiter
	always@(*) begin
		if (ao[0] == 1'b1) address = 2'd0;
		else if (ao[1] == 1'b1) address = 2'd1;
		else if (ao[2] == 1'b1) address = 2'd2;		
		//else if (ao[3] == 1'b1) address = 4'd3;
		//else if (ao[4] == 1'b1) address = 4'd4;
		//else if (ao[5] == 1'b1) address = 4'd5;
		//else if (ao[6] == 1'b1) address = 4'd6;
		//else if (ao[7] == 1'b1) address = 4'd7;
		//else if (ao[8] == 1'b1) address = 4'd8;
		//else if (ao[9] == 1'b1) address = 4'd9;
		//else if (ao[10] == 1'b1) address = 4'd10;
		//else if (ao[11] == 1'b1) address = 4'd11;
		else address = 2'd0;   // handle erroneous address encoding
	end
*/

/*
wire ao_or;
assign ao_or = ao[0] | ao[1];
// LOGIC to generate N-bit address. Assumes only one of the ao bits are asserted/ granted by the arbiter
	always@(posedge ao_or) begin
		if (ao[0] == 1'b1 && ao[1] == 1'b0) address = 1'd0;
		else if (ao[1] == 1'b1 && ao[0] == 1'b0) address = 1'd1;
		//else if (ao[2] == 1'b1) address <= 2'd2;		
		//else if (ao[3] == 1'b1) address = 4'd3;
		//else if (ao[4] == 1'b1) address = 4'd4;
		//else if (ao[5] == 1'b1) address = 4'd5;
		//else if (ao[6] == 1'b1) address = 4'd6;
		//else if (ao[7] == 1'b1) address = 4'd7;
		//else if (ao[8] == 1'b1) address = 4'd8;
		//else if (ao[9] == 1'b1) address = 4'd9;
		//else if (ao[10] == 1'b1) address = 4'd10;
		//else if (ao[11] == 1'b1) address = 4'd11;
		//else address = address;   // handle erroneous address encoding
	end
*/
	always@(*) begin
		if (ao[0] == 1'b1) address = 3'd0;
		else if (ao[1] == 1'b1) address = 3'd1;
		else if (ao[2] == 1'b1) address = 3'd2;		
		else if (ao[3] == 1'b1) address = 3'd3;
		else if (ao[4] == 1'b1) address = 3'd4;
		else if (ao[5] == 1'b1) address = 3'd5;
		else if (ao[6] == 1'b1) address = 3'd6;
		else if (ao[7] == 1'b1) address = 3'd7;
		//else if (ao[8] == 1'b1) address = 4'd8;
		//else if (ao[9] == 1'b1) address = 4'd9;
		//else if (ao[10] == 1'b1) address = 4'd10;
		//else if (ao[11] == 1'b1) address = 4'd11;
		else address = 3'd0;   // handle erroneous address encoding
	end
endmodule
