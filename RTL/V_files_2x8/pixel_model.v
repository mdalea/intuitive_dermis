`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    12:23:30 10/27/2021 
// Design Name: 
// Module Name:    pixel_model 
// Project Name: 
// Target Devices: 
// Tool versions: 
// Description: 
//
// Dependencies: 
//
// Revision: 
// Revision 0.01 - File Created
// Additional Comments: 
//
//////////////////////////////////////////////////////////////////////////////////
module pixel_model(
	input rst,

	input spk_in,

	input acky,
	input ackx_pulse,

	//output reg q_latch
	output wire q_latch
  );

	wire set;
	wire reset;
	
	assign set =  ~rst & spk_in & ~acky;
   assign reset = rst | (acky & ackx_pulse);




/*	// model gated SR-latch for spike input
	always@(*) begin
		if (rst==1) begin
			q_latch <= 1'b0;
		end else begin
			if (set == 1'b1) begin  // once row is selected, further spiking from same row is prohibited
					q_latch <= 1'b1;
			end else if (reset==1'b1) begin
					q_latch <= 1'b0;			
			end else
					q_latch <= q_latch;
		end
	end*/


sr_latch_nor u0 (.set( set ),
			.reset( reset ),
			.q( q_latch ) );


endmodule
