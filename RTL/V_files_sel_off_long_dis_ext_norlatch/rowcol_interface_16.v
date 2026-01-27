//Verilog HDL for "spike_enc", "rowcol_interface" "functional"


module rowcol_interface_16 (

input rst,
input greedy,
input [15:0] n_p,
input n_ai,
input [15:0] n_ri,

input arbtop_n_ri,  // added so spike is not discarded

output wire [15:0] ro,
output wire [15:0] s,
output wire [15:0] ao,
input aer_dis 
);


rowcol_interface r0 ( 
			.rst( rst ), .greedy(greedy),
			.n_p( n_p[0] ),
			.n_ai ( n_ai ),
			.n_ri ( n_ri[0] ),
			
			.arbtop_n_ri(arbtop_n_ri),			
			
			.ro ( ro[0] ),
			.s ( s[0] ),
			.ao ( ao[0] )
			, .aer_dis ( aer_dis )     );
			
rowcol_interface r1( 
			.rst( rst ), .greedy(greedy),
			.n_p( n_p[1] ),
			.n_ai ( n_ai ),
			.n_ri ( n_ri[1] ),
			
			.arbtop_n_ri(arbtop_n_ri),				
			
			.ro ( ro[1] ),
			.s ( s[1] ),
			.ao ( ao[1] )
			, .aer_dis ( aer_dis )     );			
			
rowcol_interface r2( 
			.rst( rst ), .greedy(greedy),
			.n_p( n_p[2] ),
			.n_ai ( n_ai ),
			.n_ri ( n_ri[2] ),

			.arbtop_n_ri(arbtop_n_ri),	
			
			.ro ( ro[2] ),
			.s ( s[2] ),
			.ao ( ao[2] )
			, .aer_dis ( aer_dis )     );			

rowcol_interface r3( 
			.rst( rst ), .greedy(greedy),
			.n_p( n_p[3] ),
			.n_ai ( n_ai ),
			.n_ri ( n_ri[3] ),

			.arbtop_n_ri(arbtop_n_ri),	
			
			.ro ( ro[3] ),
			.s ( s[3] ),
			.ao ( ao[3] )
			, .aer_dis ( aer_dis )     );		

rowcol_interface r4( 
			.rst( rst ), .greedy(greedy),
			.n_p( n_p[4] ),
			.n_ai ( n_ai ),
			.n_ri ( n_ri[4] ),

			.arbtop_n_ri(arbtop_n_ri),	
			
			.ro ( ro[4] ),
			.s ( s[4] ),
			.ao ( ao[4] )
			, .aer_dis ( aer_dis )     );					

rowcol_interface r5( 
			.rst( rst ), .greedy(greedy),
			.n_p( n_p[5] ),
			.n_ai ( n_ai ),
			.n_ri ( n_ri[5] ),

			.arbtop_n_ri(arbtop_n_ri),	
			
			.ro ( ro[5] ),
			.s ( s[5] ),
			.ao ( ao[5] )
			, .aer_dis ( aer_dis )     );								

rowcol_interface r6( 
			.rst( rst ), .greedy(greedy),
			.n_p( n_p[6] ),
			.n_ai ( n_ai ),
			.n_ri ( n_ri[6] ),

			.arbtop_n_ri(arbtop_n_ri),	
			
			.ro ( ro[6] ),
			.s ( s[6] ),
			.ao ( ao[6] )
			, .aer_dis ( aer_dis )     );								
			

rowcol_interface r7( 
			.rst( rst ), .greedy(greedy),
			.n_p( n_p[7] ),
			.n_ai ( n_ai ),
			.n_ri ( n_ri[7] ),

			.arbtop_n_ri(arbtop_n_ri),	
			
			.ro ( ro[7] ),
			.s ( s[7] ),
			.ao ( ao[7] )
			, .aer_dis ( aer_dis )     );				


rowcol_interface r8( 
			.rst( rst ), .greedy(greedy),
			.n_p( n_p[8] ),
			.n_ai ( n_ai ),
			.n_ri ( n_ri[8] ),

			.arbtop_n_ri(arbtop_n_ri),	
			
			.ro ( ro[8] ),
			.s ( s[8] ),
			.ao ( ao[8] )
			, .aer_dis ( aer_dis )     );				
			

rowcol_interface r9( 
			.rst( rst ), .greedy(greedy),
			.n_p( n_p[9] ),
			.n_ai ( n_ai ),
			.n_ri ( n_ri[9] ),

			.arbtop_n_ri(arbtop_n_ri),	
			
			.ro ( ro[9] ),
			.s ( s[9] ),
			.ao ( ao[9] )
			, .aer_dis ( aer_dis )     );

rowcol_interface r10( 
	.rst( rst ), .greedy(greedy),
	.n_p( n_p[10] ),
	.n_ai ( n_ai ),
	.n_ri ( n_ri[10] ),

	.arbtop_n_ri(arbtop_n_ri),

	.ro ( ro[10] ),
	.s ( s[10] ),
	.ao ( ao[10] )
, .aer_dis ( aer_dis )     );

rowcol_interface r11( 
	.rst( rst ), .greedy(greedy),
	.n_p( n_p[11] ),
	.n_ai ( n_ai ),
	.n_ri ( n_ri[11] ),

	.arbtop_n_ri(arbtop_n_ri),

	.ro ( ro[11] ),
	.s ( s[11] ),
	.ao ( ao[11] )
, .aer_dis ( aer_dis )     );

rowcol_interface r12( 
	.rst( rst ), .greedy(greedy),
	.n_p( n_p[12] ),
	.n_ai ( n_ai ),
	.n_ri ( n_ri[12] ),

	.arbtop_n_ri(arbtop_n_ri),

	.ro ( ro[12] ),
	.s ( s[12] ),
	.ao ( ao[12] )
, .aer_dis ( aer_dis )     );

rowcol_interface r13( 
	.rst( rst ), .greedy(greedy),
	.n_p( n_p[13] ),
	.n_ai ( n_ai ),
	.n_ri ( n_ri[13] ),

	.arbtop_n_ri(arbtop_n_ri),

	.ro ( ro[13] ),
	.s ( s[13] ),
	.ao ( ao[13] )
, .aer_dis ( aer_dis )     );

rowcol_interface r14( 
	.rst( rst ), .greedy(greedy),
	.n_p( n_p[14] ),
	.n_ai ( n_ai ),
	.n_ri ( n_ri[14] ),

	.arbtop_n_ri(arbtop_n_ri),

	.ro ( ro[14] ),
	.s ( s[14] ),
	.ao ( ao[14] )
, .aer_dis ( aer_dis )     );

rowcol_interface r15( 
	.rst( rst ), .greedy(greedy),
	.n_p( n_p[15] ),
	.n_ai ( n_ai ),
	.n_ri ( n_ri[15] ),

	.arbtop_n_ri(arbtop_n_ri),

	.ro ( ro[15] ),
	.s ( s[15] ),
	.ao ( ao[15] )
, .aer_dis ( aer_dis )     );

			
						
endmodule
