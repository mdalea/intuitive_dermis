`timescale 1ns / 1ps

////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer:
//
// Create Date:   19:03:15 10/26/2021
// Design Name:   aer_top
// Module Name:   /users/micas/malea/RTL/aer_interface/tb_aer_top.v
// Project Name:  aer_interface
// Target Device:  
// Tool versions:  
// Description: 
//
// Verilog Test Fixture created by ISE for module: aer_top
//
// Dependencies:
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
////////////////////////////////////////////////////////////////////////////////

module tb_aer_top;

	// Inputs to AER
	reg ack;
	reg greedy;
	wire [11:0] n_p;
	wire [15:0] n_cox_on;
	wire [15:0] n_cox_off;
	wire [15:0] cix_on_del;
	wire [15:0] cix_off_del;	
	
	reg aer_dis;

	// Outputs to AER
	wire reqon;
	wire reqoff;
	wire [3:0] addrx;
	wire [3:0] addry;
	wire [11:0] s;
	wire [15:0] cix_on;
	wire [15:0] cix_off;
	
	reg [10:0] s_ext;
	reg rst;
	reg [191:0] spk_in_on;
	reg [191:0] spk_in_off;

	// Instantiate the Unit Under Test (UUT)
	aer_row_off_long_dis_ext_norlatch row_uut (
		.RST(rst),
		.GREEDY(greedy),
		.ACK(ack), 
		.ADDRY(addry), 
		.N_P(n_p),  
		.S(s),
		.AER_DIS(aer_dis),
		
		.S_EXT(s_ext)
	);

	aer_col_off_long_dis_norlatch col_uut (
		.RST(rst),
		.GREEDY(greedy),
		.ACK(ack), 
		.REQON(reqon), 
		.REQOFF(reqoff), 
		.ADDRX(addrx), 
		.N_COX_ON(n_cox_on), 
		.N_COX_OFF(n_cox_off), 		
		.CIX_ON(cix_on), 
		.CIX_OFF(cix_off), 
		.CIX_ON_DEL(cix_on_del),
		.CIX_OFF_DEL(cix_off_del),
		.AER_DIS(aer_dis)
	);

	pixel_arr_model pix_arr(
		.rst(rst),
		.spk_in_on(spk_in_on),
		.spk_in_off(spk_in_off),

		.s(s),
		.cix_on(cix_on),
		.cix_off(cix_off),
		.n_p(n_p),
		.n_cox_on(n_cox_on),
		.n_cox_off(n_cox_off),

		.cix_on_del(cix_on_del),		
		.cix_off_del(cix_off_del)
    );
	 
	 always@(*) begin
		//if (reqon==1 || reqoff==1) begin
		//	ack <= # 1;
		//end else if (reqon==0 || reqoff==0) begin
		////	ack <= # 0;
		//end
		if (rst==1) 
			ack <= 0;
		else
			ack <= #10 (reqon | reqoff); 
	 end
	 
	initial begin
		// Initialize Inputs
		
		
		aer_dis=0;
		s_ext=0;
		
		rst=0;
		#100;
		//spk_in = 100'b0;
		//ack=1;  // serves as the reset by resetting the rowcol interface REQs to 0, propagates back to ROWSEL (s), enabling gated SR latch
		//#1000;
		//ack=0;   
		rst=1;
		#1000;
		rst=0;
		#100;
		rst=1;
		#1000;
		rst=0;
		greedy =0;
		
		#100;
		spk_in_on[6] = 1;
		spk_in_on[66] = 1;
		#1000;
		spk_in_on[6] = 0;
		spk_in_on[66] = 0;
		#100;
		


		spk_in_on[191:0] = ~0;
		#1000;
		spk_in_on[191:0] = 0;
		
		#100000;
		spk_in_off[191:0] = ~0;
		#1000;
		spk_in_off[191:0] = 0;
		
/*		spk_in[99:90] = ~0; // will not be serviced until ack for previous REQ is received
		#10
		spk_in[99:90] = 0;
		#8000;
		spk_in[9:0] = ~0;
		#10;
		spk_in[9:0] = 0;*/

		

		#1000;

		spk_in_on[0] = ~0; // will not be serviced until ack for previous REQ is received
		spk_in_on[69] = ~0;
		spk_in_on[49] = ~0;		
		spk_in_on[191] = ~0;
		
//		//spk_in_off[0] = ~0; // will not be serviced until ack for previous REQ is received	
//		//spk_in_off[69] = ~0;
//		//spk_in_off[49] = ~0;		
//		//spk_in_off[191] = ~0;	
//
		#100000;
//
//		//spk_in_on[0] = ~0; // will not be serviced until ack for previous REQ is received
//		//spk_in_on[69] = ~0;
//		//spk_in_on[49] = ~0;		
//		//spk_in_on[191] = ~0;
//		
		spk_in_off[0] = ~0; // will not be serviced until ack for previous REQ is received	
		spk_in_off[69] = ~0;
		spk_in_off[49] = ~0;		
		spk_in_off[191] = ~0;			

		#50000;
		spk_in_on[0] = 0;
		spk_in_on[69] = 0;
		spk_in_on[49] = 0;			
		spk_in_on[191] = 0;
		
		//spk_in_off[0] = 0; // will not be serviced until ack for previous REQ is received		
		//spk_in_off[69] = 0;
		//spk_in_off[49] = 0;			
		//spk_in_off[191] = 0;
		
		#50000;
//		//spk_in_on[0] = 0;
//		//spk_in_on[69] = 0;
//		//spk_in_on[49] = 0;			
//		//spk_in_on[191] = 0;
//		
		spk_in_off[0] = 0; // will not be serviced until ack for previous REQ is received		
		spk_in_off[69] = 0;
		spk_in_off[49] = 0;			
		spk_in_off[191] = 0;
//		
//		#100000;
//		greedy=1;
//		spk_in_on[0] = ~0; // will not be serviced until ack for previous REQ is received
//		spk_in_on[191] = ~0;
//		#50000;
//		spk_in_on[0] = 0;
//		spk_in_on[191] = 0;     

//    ON/OFF spikes from different address

		#100000;

		spk_in_on[0] = ~0; // will not be serviced until ack for previous REQ is received
		spk_in_on[69] = ~0;
		spk_in_on[49] = ~0;		
		spk_in_on[191] = ~0;
		
//		//spk_in_off[0] = ~0; // will not be serviced until ack for previous REQ is received	
//		//spk_in_off[69] = ~0;
//		//spk_in_off[49] = ~0;		
//		//spk_in_off[191] = ~0;	
//
		#100000;
//
//		//spk_in_on[0] = ~0; // will not be serviced until ack for previous REQ is received
//		//spk_in_on[69] = ~0;
//		//spk_in_on[49] = ~0;		
//		//spk_in_on[191] = ~0;
//		
		spk_in_off[10] = ~0; // will not be serviced until ack for previous REQ is received	
		spk_in_off[79] = ~0;
		spk_in_off[59] = ~0;		
		spk_in_off[181] = ~0;			

		#50000;
		spk_in_on[0] = 0;
		spk_in_on[69] = 0;
		spk_in_on[49] = 0;			
		spk_in_on[191] = 0;
		
		//spk_in_off[0] = 0; // will not be serviced until ack for previous REQ is received		
		//spk_in_off[69] = 0;
		//spk_in_off[49] = 0;			
		//spk_in_off[191] = 0;
		
		#50000;
//		//spk_in_on[0] = 0;
//		//spk_in_on[69] = 0;
//		//spk_in_on[49] = 0;			
//		//spk_in_on[191] = 0;
//		
		spk_in_off[10] = 0; // will not be serviced until ack for previous REQ is received		
		spk_in_off[79] = 0;
		spk_in_off[59] = 0;			
		spk_in_off[181] = 0;
//		
//		#100000;
//		greedy=1;
//		spk_in_on[0] = ~0; // will not be serviced until ack for previous REQ is received
//		spk_in_on[191] = ~0;
//		#50000;
//		spk_in_on[0] = 0;
//		spk_in_on[191] = 0;  



/*******************/
		aer_dis=1;
		spk_in_on[191:0] = ~0;
		#1000;
		spk_in_on[191:0] = 0;
		
		#100000;
		spk_in_off[191:0] = ~0;
		#1000;
		spk_in_off[191:0] = 0;
		
/*		spk_in[99:90] = ~0; // will not be serviced until ack for previous REQ is received
		#10
		spk_in[99:90] = 0;
		#8000;
		spk_in[9:0] = ~0;
		#10;
		spk_in[9:0] = 0;*/

		

		#1000;

		spk_in_on[0] = ~0; // will not be serviced until ack for previous REQ is received
		spk_in_on[69] = ~0;
		spk_in_on[49] = ~0;		
		spk_in_on[191] = ~0;
		
//		//spk_in_off[0] = ~0; // will not be serviced until ack for previous REQ is received	
//		//spk_in_off[69] = ~0;
//		//spk_in_off[49] = ~0;		
//		//spk_in_off[191] = ~0;	
//
		#100000;
//
//		//spk_in_on[0] = ~0; // will not be serviced until ack for previous REQ is received
//		//spk_in_on[69] = ~0;
//		//spk_in_on[49] = ~0;		
//		//spk_in_on[191] = ~0;
//		
		spk_in_off[0] = ~0; // will not be serviced until ack for previous REQ is received	
		spk_in_off[69] = ~0;
		spk_in_off[49] = ~0;		
		spk_in_off[191] = ~0;			

		#50000;
		spk_in_on[0] = 0;
		spk_in_on[69] = 0;
		spk_in_on[49] = 0;			
		spk_in_on[191] = 0;
		
		//spk_in_off[0] = 0; // will not be serviced until ack for previous REQ is received		
		//spk_in_off[69] = 0;
		//spk_in_off[49] = 0;			
		//spk_in_off[191] = 0;
		
		#50000;
//		//spk_in_on[0] = 0;
//		//spk_in_on[69] = 0;
//		//spk_in_on[49] = 0;			
//		//spk_in_on[191] = 0;
//		
		spk_in_off[0] = 0; // will not be serviced until ack for previous REQ is received		
		spk_in_off[69] = 0;
		spk_in_off[49] = 0;			
		spk_in_off[191] = 0;
//		
//		#100000;
//		greedy=1;
//		spk_in_on[0] = ~0; // will not be serviced until ack for previous REQ is received
//		spk_in_on[191] = ~0;
//		#50000;
//		spk_in_on[0] = 0;
//		spk_in_on[191] = 0;     

//    ON/OFF spikes from different address

		#100000;

		spk_in_on[0] = ~0; // will not be serviced until ack for previous REQ is received
		spk_in_on[69] = ~0;
		spk_in_on[49] = ~0;		
		spk_in_on[191] = ~0;
		
//		//spk_in_off[0] = ~0; // will not be serviced until ack for previous REQ is received	
//		//spk_in_off[69] = ~0;
//		//spk_in_off[49] = ~0;		
//		//spk_in_off[191] = ~0;	
//
		#100000;
//
//		//spk_in_on[0] = ~0; // will not be serviced until ack for previous REQ is received
//		//spk_in_on[69] = ~0;
//		//spk_in_on[49] = ~0;		
//		//spk_in_on[191] = ~0;
//		
		spk_in_off[10] = ~0; // will not be serviced until ack for previous REQ is received	
		spk_in_off[79] = ~0;
		spk_in_off[59] = ~0;		
		spk_in_off[181] = ~0;			

		#50000;
		spk_in_on[0] = 0;
		spk_in_on[69] = 0;
		spk_in_on[49] = 0;			
		spk_in_on[191] = 0;
		
		//spk_in_off[0] = 0; // will not be serviced until ack for previous REQ is received		
		//spk_in_off[69] = 0;
		//spk_in_off[49] = 0;			
		//spk_in_off[191] = 0;
		
		#50000;
//		//spk_in_on[0] = 0;
//		//spk_in_on[69] = 0;
//		//spk_in_on[49] = 0;			
//		//spk_in_on[191] = 0;
//		
		spk_in_off[10] = 0; // will not be serviced until ack for previous REQ is received		
		spk_in_off[79] = 0;
		spk_in_off[59] = 0;			
		spk_in_off[181] = 0;
//		
//		#100000;
//		greedy=1;
//		spk_in_on[0] = ~0; // will not be serviced until ack for previous REQ is received
//		spk_in_on[191] = ~0;
//		#50000;
//		spk_in_on[0] = 0;
//		spk_in_on[191] = 0;  

   
		// Add stimulus here	
		// Add stimulus here

	end
      
endmodule

