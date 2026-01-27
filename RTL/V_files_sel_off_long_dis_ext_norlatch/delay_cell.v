// model generation of both spkdel and spkpulse

`timescale 1ns / 1ps


module delay_cell(
	input cix,
	output reg cix_del
    );



	// model cix (ACKX) delay for spike height generation
	always@(*) begin
		//if (cix == 1)   #10000 cix_del =  cix;
		if (cix == 1)    cix_del = #500 cix;
		else if (cix == 0)  cix_del =  cix;
	end



endmodule
