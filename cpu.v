

///////////////////////////////// if/id pipeline register ///////////////////////////////////////////
module if_id_reg (inst,pc,inst_out,pc_out,taken,hold,jump,jal,jr,clk);

input [31:0] inst;
input [31:0] pc;
input taken,hold,jump,jal,jr;
input clk;
output [31:0] inst_out;
output [31:0] pc_out;

reg [0:31] if_id_cont [0:1];  //0:inst  1:pc_4

always@(posedge clk)
begin
/* if the branch is taken so the instruction in this stage should be nop also for jump but hold means we have to keep the same instruction*/
if(taken==1 | jump==1 | jal==1 | jr==1)
	begin
		if_id_cont[0]<=32'b0;
		if_id_cont[1]<=pc;
	end
else if(hold==1)
	begin
		if_id_cont[0]<=if_id_cont[0];
		if_id_cont[1]<=if_id_cont[1];
	end
else
	begin
		if_id_cont[0]<=inst;
		if_id_cont[1]<=pc;
	end
end

assign inst_out=if_id_cont[0];
assign pc_out=if_id_cont[1];

initial
begin
	if_id_cont[0]=32'b0;
	if_id_cont[1]=32'b0;
end
endmodule

//////////////////////// id/ex pipeline register  /////////////////////////////////////////
module id_ex_reg (ctrl_signals,pc,read_data1,read_data2,after_sign_ext,rd,rt,rs,
			ctrl_signals_out,pc_out,read_data1_out,read_data2_out,after_sign_ext_out,rd_out,rt_out,rs_out,flush,taken,clk);

input [10:0] ctrl_signals;
input [31:0] pc;
input [31:0] read_data1;
input [31:0] read_data2;
input [31:0] after_sign_ext;
input [4:0] rd;
input [4:0] rt;
input [4:0] rs;
input flush,taken;
input clk;

output [10:0] ctrl_signals_out;
output [31:0] pc_out;
output [31:0] read_data1_out;
output [31:0] read_data2_out;
output [31:0] after_sign_ext_out;
output [4:0] rd_out;
output [4:0] rt_out;
output [4:0] rs_out;

reg [10:0] ctrl_signals_reg;
reg [31:0] data_reg[0:3]; // 0:read_data1   1:read_data2  2:pc  3:after_sign_ext
reg [4:0] reg_regs [0:2]; // 0:rt  1:rd  2: rs

initial
begin
ctrl_signals_reg=11'b0;
data_reg[0]=32'b0;data_reg[1]=32'b0;data_reg[2]=32'b0;data_reg[3]=32'b0;
reg_regs[0]=5'b0;reg_regs[1]=5'b0;reg_regs[2]=5'b0;
end

/*here branch is determined so if taken the ctrl_signals will be zero so for flush means that the instruction will give wrong outputs*/
always@(posedge clk)
begin
if(flush ==1 | taken==1) ctrl_signals_reg<=11'b0;
else ctrl_signals_reg<=ctrl_signals;

data_reg[0]<=read_data1;
data_reg[1]<=read_data2;
data_reg[2]<=pc;
data_reg[3]<=after_sign_ext;

reg_regs[0]<=rt;
reg_regs[1]<=rd;
reg_regs[2]<=rs;

end

assign ctrl_signals_out=ctrl_signals_reg;

assign read_data1_out=data_reg[0];
assign read_data2_out=data_reg[1];
assign pc_out=data_reg[2];
assign after_sign_ext_out=data_reg[3];

assign rt_out=reg_regs[0];
assign rd_out=reg_regs[1];
assign rs_out=reg_regs[2];

endmodule

///////////////////////////////// ex/mem pipeline register /////////////////////////////////
module ex_mem_reg (ctrl_signals,alu_out,read_data2,reg_dst_out,
			ctrl_signals_out,alu_out_out,read_data2_out,reg_dst_out_out,clk);

input [3:0] ctrl_signals;
input [4:0] reg_dst_out;
input [31:0] alu_out;
input [31:0] read_data2;
input clk;

output [3:0] ctrl_signals_out;
output [4:0] reg_dst_out_out;
output [31:0] alu_out_out;
output [31:0] read_data2_out;

reg [3:0] ctrl_signals_reg;
reg [4:0] reg_dst_reg;
reg [31:0] data_regs[0:1]; // 0:alu_out  1:read_data2

initial
begin
ctrl_signals_reg=4'b0;
reg_dst_reg=5'b0;
data_regs[0]=32'b0;
data_regs[1]=32'b0;
end

always@(posedge clk)
begin
ctrl_signals_reg<=ctrl_signals;
reg_dst_reg<=reg_dst_out;
data_regs[0]<=alu_out;
data_regs[1]<=read_data2;
end

assign ctrl_signals_out=ctrl_signals_reg;
assign reg_dst_out_out=reg_dst_reg;
assign alu_out_out=data_regs[0];
assign read_data2_out=data_regs[1];

endmodule

/////////////////////////// mem/wb pipeline register  ///////////////////////////
module mem_wb_reg(ctrl_signals,reg_dst_out,read_data,alu_out,
			ctrl_signals_out,reg_dst_out_out,read_data_out,alu_out_out,clk);

input [1:0] ctrl_signals;
input [4:0] reg_dst_out;
input [31:0] read_data;
input [31:0] alu_out;
input clk;

output [1:0] ctrl_signals_out;
output [4:0] reg_dst_out_out;
output [31:0] read_data_out;
output [31:0] alu_out_out;

reg [1:0] ctrl_signals_reg;
reg [4:0] reg_dst_reg;
reg [31:0] data_regs[0:1]; // 0: read_data   1: alu_out

initial
begin
ctrl_signals_reg=2'b0;
reg_dst_reg=5'b0;
data_regs[0]=32'b0;
data_regs[1]=32'b0;
end

always@(posedge clk)
begin

ctrl_signals_reg<=ctrl_signals;
reg_dst_reg<=reg_dst_out;
data_regs[0]<=read_data;
data_regs[1]<=alu_out;

end

assign ctrl_signals_out=ctrl_signals_reg;
assign reg_dst_out_out=reg_dst_reg;
assign read_data_out=data_regs[0];
assign alu_out_out=data_regs[1];

endmodule

///////////////// pc_counter module to determine the next pc ////////////////////////////////
module pc_counter(pc,pc_id_exe,pc_if_id,after_sign_ext_id_exe,taken,jump,jal,jr,inst_if_id,jr_cont,hold,clk);

output reg  [31:0] pc;
input [31:0] after_sign_ext_id_exe;
input [31:0] pc_id_exe;
input [31:0] pc_if_id;
input taken,jump,jal,jr;
input [31:0] inst_if_id;
input [31:0] jr_cont;
wire [31:0] mux2out;
wire [31:0] next_pc;
input hold;
input clk;

wire [31:0] pc_4;
wire [31:0] pc_4_branch;
wire [31:0] pc_branch;
wire [31:0] pc_jump;
wire [31:0] mux1out;

adder adder1(pc,32'b100,pc_4);
adder adder2(pc_id_exe,32'b100,pc_4_branch);
adder adder3(pc_if_id,32'b100,pc_jump);
adder adder4(pc_4_branch,after_sign_ext_id_exe<<2,pc_branch);
mux2to1 mux1(pc_4,pc_branch,taken,mux1out); // if the bne or beq is true so pc will be the offset+pc+4
mux2to1 mux2(mux1out,{{pc_jump[31:28]},{inst_if_id[25:0]},2'b0},jump|jal,mux2out); // if jump or jal is true so pc will be the label
mux2to1 mux3(mux2out,jr_cont,jr,next_pc); // if jump register is true so pc will be the content of thr register

initial
begin
pc=0;
end
always@(posedge clk)
if(hold==1) pc<=pc;
else pc<=next_pc;

endmodule

////////////////////////// instruction_memory ///////////////////////
module inst_mem (read_address,instruction);
input [31:0] read_address;
output wire [31:0] instruction;
reg [31:0] inst_memory[0:1023];

initial
begin
// change the directory to where you put the files in your pc
$readmemb("E:/courses/co/project/inst_mem.txt",inst_memory);
end

assign instruction=inst_memory[read_address>>2];

endmodule

/////////////////////// data_memory ///////////////////////
module data_mem (address, write_data , memRead ,memWrite,read_data,clk);
input [31:0] address;
input [31:0] write_data ;
input memRead ;
input memWrite;
input clk;
output wire [31:0] read_data;
reg [31:0] data_memory[0:1023];

initial
begin
// change the directory to where you put the files in your pc
$readmemb("E:/courses/co/project/data_mem.txt",data_memory);
end

always@(posedge clk)
begin
if(memWrite)
	data_memory[address>>2]<=write_data;
end

assign read_data=data_memory[address>>2];

endmodule

//////////////////////// register file /////////////////////////////
module reg_file (read_reg1,read_reg2,write_reg,write_data,reg_write,read_data1,read_data2,jal,pcplus4,clk);
input [4:0] read_reg1;
input [4:0] read_reg2;
input [4:0] write_reg;
input [31:0] write_data;
input reg_write;
input jal;
input [31:0]pcplus4;
input clk;
output wire [31:0] read_data1;
output wire [31:0] read_data2;
reg [31:0] reg_f[0:31];
initial
begin
// change the directory to where you put the files in your pc
$readmemb("E:/courses/co/project/reg_f.txt",reg_f);
end

assign read_data1=reg_f[read_reg1];
assign read_data2=reg_f[read_reg2];

always@(negedge clk)
begin
if(reg_write)
	reg_f[write_reg]<=write_data;
	
if (jal && write_reg!=31) // jal is determined in the decoding stage so it writes the value of pc+4 in ra register
   	reg_f[31]<=pcplus4;
end

endmodule

///////////////////////// alu ////////////////////////////////////////////
module real_alu(in1,in2,shamt,op,out,zero_flag);
input [31:0] in1;
input [31:0] in2;
input [3:0] op;
input [4:0] shamt;
output reg  [31:0] out;
output reg zero_flag;
parameter And=4'b0000/*0*/, Or=4'b0001/*1*/, Add=4'b0010/*2*/, Sub=4'b0110/*6*/, Slt=4'b0111/*7*/, Nor=4'b1100/*12*/, Sll=4'b1101 /*13*/,
	  Srl=4'b0011/*3*/, Sra=4'b0100/*4*/, lui=4'b1010/*10*/, slti=4'b0101;
always@(in1,in2,op,shamt)
begin
zero_flag<=0;
case(op)
	And: out<=(in1&in2);
	Or: out<=(in1|in2);
	Add: out<=(in1+in2);
	Sub:
		begin 
		out<=(in1-in2);
		if(in1==in2) zero_flag<=1;
		else zero_flag<=0;
		end

	//slt: the value of out=the sign bit of in1-in2
	
	Slt:
		begin
		out<=(32'b1 & ((in1-in2)>>31));
		end

	Nor: out<=(~(in1|in2));
	Sll: out<=(in2<<shamt);
	Srl: out<=(in2>>shamt);
	Sra: out<=($signed(in2)>>>shamt);
	lui: out<=(in2<<16);
	slti: 
		begin
		out<=(32'b1 & ((in1-in2)>>31));
		end
	
	default:  out<=32'b0;

endcase
end

endmodule

/////////////////////// hazard detection unit ////////////////////////////
module HDU(id_ex_mem_read,id_ex_rt,if_id_rs,if_id_rt,hold);

input id_ex_mem_read;
input [4:0] id_ex_rt;
input [4:0] if_id_rs;
input [4:0] if_id_rt;

output reg hold;

always@(id_ex_mem_read,id_ex_rt,if_id_rs,if_id_rt)
begin
if(id_ex_mem_read & ((id_ex_rt==if_id_rs) | (id_ex_rt==if_id_rt)))  hold=1;
else  hold=0;
end
endmodule

////////////////////////// control unit ///////////////////////////
module ctrl_unit (opcode,mem_write,reg_dst,jump ,branch,mem_to_reg,mem_read,alu_src,reg_write,aluop,
branch_not_eq,jal, jr,  function_bits);
input [5:0] opcode;
input [5:0] function_bits;  // to determine jr in the decoding stage
output reg mem_write;
output reg reg_dst;
output reg jump;
output reg branch;
output reg branch_not_eq; //added to support bne
output reg mem_to_reg;
output reg mem_read;
output reg alu_src;
output reg reg_write;
output reg [2:0] aluop; //aluop became 3 bits not 2 to support more operations
output reg jal;
output reg jr;

always@(opcode,function_bits)
begin
case(opcode)

//R opcode=0
	6'b000_000:
         case(function_bits)

			//jump register: opcode=0; fn=8; it jumps to the address carried by $ra
			//it wont write in RF //same signals as jump except jump=0

		6'b001000: {mem_write,reg_dst,jump ,branch,branch_not_eq,mem_to_reg,mem_read,alu_src,reg_write,aluop,  jal, jr}	<=14'b0x000x0x0xxx01;

			//other R instructions

		default: {mem_write,reg_dst,jump ,branch,branch_not_eq,mem_to_reg,mem_read,alu_src,reg_write,aluop,  jal, jr}	<=14'b01000000101000;

	endcase


//jump: opcode=2

	6'b000010: {mem_write,reg_dst,jump ,branch,branch_not_eq,mem_to_reg,mem_read,alu_src,reg_write,aluop,  jal, jr}<=14'b0x100x0x0xxx00;

//jump and link: opcode=3, jumps and saves address of following instr in $ra
//jump signal will be zero, jal will be one, RF will write=1

	6'b000011: {mem_write,reg_dst,jump ,branch,branch_not_eq,mem_to_reg,mem_read,alu_src,reg_write,aluop,  jal, jr}<=14'b0x000x0x1xxx10;


//lw: opcode=35

	6'b100011: {mem_write,reg_dst,jump ,branch,branch_not_eq,mem_to_reg,mem_read,alu_src,reg_write,aluop,  jal, jr}<=14'b00000111100000;

//sw: opcode=43

	6'b101011: {mem_write,reg_dst,jump ,branch,branch_not_eq,mem_to_reg,mem_read,alu_src,reg_write,aluop,  jal, jr}<=14'b1x000x01000000;

//beq: opcode=4

	6'b000100: {mem_write,reg_dst,jump ,branch,branch_not_eq,mem_to_reg,mem_read,alu_src,reg_write,aluop,  jal, jr}<=14'b0x010x00000100;

//bne:op code=5
//same bits as beq except for branch will be zero and branch not will be one

	6'b000101: {mem_write,reg_dst,jump ,branch,branch_not_eq,mem_to_reg,mem_read,alu_src,reg_write,aluop,  jal, jr}<=14'b0x001x00000100;
	

//addi: op=8
//doesnot write in memory, writes in rt, no jump or beq or bne, mem to reg=0 to take alu o/p, doesnot read from mem
//alu src=1 to take immediate value, writes in RF=1, aluop=000 to add

	6'b001000: {mem_write,reg_dst,jump ,branch,branch_not_eq,mem_to_reg,mem_read,alu_src,reg_write,aluop,  jal, jr}<=14'b00000001100000;

//ori: op=13 :1101
//same bits as addi except in alu_op=100 for or operation

	6'b001101: {mem_write,reg_dst,jump ,branch,branch_not_eq,mem_to_reg,mem_read,alu_src,reg_write,aluop,  jal, jr}<=14'b00000001110000;

//andi:op=12: 1100
//same bits as ori except in alu_op=011

	6'b001100: {mem_write,reg_dst,jump ,branch,branch_not_eq,mem_to_reg,mem_read,alu_src,reg_write,aluop,  jal, jr}<=14'b00000001101100;

//slti: op=10 : 1010
//same bits as ori but alu will slti, alu_op=110

	6'b001010: {mem_write,reg_dst,jump ,branch,branch_not_eq,mem_to_reg,mem_read,alu_src,reg_write,aluop,  jal, jr}<=14'b00000001111000;

//lui: op=15: 1111
//takes the immediate value then shifts it by 16 to left, stores it in $rt
//doesnot write ine mem, writes in rt,no jumo, no beq.no bne, memtoreg=0 it takes value of alu o/p,
//doesnot read from mem, aku_src=1 uses immediate value, writes in reg, uses alu_op=101 to shift

	6'b001111: {mem_write,reg_dst,jump ,branch,branch_not_eq,mem_to_reg,mem_read,alu_src,reg_write,aluop, jal, jr}<=14'b00000001110100;

	default: {mem_write,reg_dst,jump ,branch,branch_not_eq,mem_to_reg,mem_read,alu_src,reg_write,aluop, jal, jr}<=14'b00000000000000;


endcase
end

endmodule

///////////////////// 32-bits 2to1 mux //////////////////////
module mux2to1(in1,in2,sel,out);
input [31:0] in1;
input [31:0] in2;
input sel;

output reg [31:0] out;

always@(in1,in2,sel)
begin
case(sel)
	1'b0: out<=in1;
	1'b1: out<=in2;
endcase
end

endmodule

//////////////////// 5-bits 2to1 mux //////////////////////////
module mux5_2to1(in1,in2,sel,out);
input [4:0] in1;
input [4:0] in2;
input sel;

output reg [4:0] out;

always@(in1,in2,sel)
begin
case(sel)
	1'b0: out<=in1;
	1'b1: out<=in2;
endcase
end

endmodule
//////////////////////////// 32-bits input 3to1 mux /////////////////////////
module mux3to1(in1,in2,in3,sel,out);
input [31:0] in1;
input [31:0] in2;
input [31:0] in3;
input [1:0] sel;
output reg [31:0] out;
always@(in1,in2,in3,sel)
begin
case(sel)
	2'b00: out<=in1;
	2'b01: out<=in2;
	2'b10: out<=in3;
	2'b11: out<=in1;
endcase
end

endmodule

////////////////////////////// adder ////////////////////////////////
module adder(in1,in2,out);
input [31:0] in1;
input [31:0] in2;

output reg [31:0] out;

always@(in1,in2)
out<=in1+in2;

endmodule

///////////////////// sign extension ///////////////////////////////
module sign_ext(in,out);
input [15:0] in;
output reg [31:0] out;
always@(in)
out<={{16{in[15]}},in}; 
endmodule

/////////////////////// alu control /////////////////////////// 
module alu_control(control_bits, aluop, function_bits );

output reg [3:0]control_bits; //extended alu_op to 3 bits instead of 2 to accomodate ori, andi

input [2:0]aluop;
input [5:0]function_bits;
parameter And=4'b0000/*0*/, Or=4'b0001/*1*/, Add=4'b0010/*2*/, Sub=4'b0110/*6*/, Slt=4'b0111/*7*/, Nor=4'b1100/*12*/, Sll=4'b1101 /*13*/,
	  Srl=4'b0011/*3*/, Sra=4'b0100/*4*/, lui=4'b1010/*10*/, slti=4'b0101; //these shifts control bits aren't standardized

always@(aluop, function_bits)
begin
case(aluop)
	3'b000: control_bits<=Add;
	3'b001: control_bits<=Sub;
	3'b010:
	begin
	case(function_bits)
		6'b100_000: control_bits<=Add;
		6'b100_010: control_bits<=Sub;
		6'b100_100: control_bits<=And;
		6'b100_101: control_bits<=Or;
		6'b101_010: control_bits<=Slt;
		6'b100_111: control_bits<=Nor;
		6'b000_000: control_bits<=Sll;
		6'b000_010: control_bits<=Srl;
		6'b000_011: control_bits<=Sra;
		default   : control_bits<=4'b0000;
	endcase
	end
	
	3'b011: control_bits<=And; //andi operation
	3'b100: control_bits<=Or;  //ori operation
	3'b101: control_bits<=lui; //load upper immediate
	3'b110: control_bits<=slti; //set on less than immediate
	default   : control_bits<=4'b0000;
endcase

end
endmodule

///////////////Forwarding unit detection module/////////////////
module forwarding_unit (forwarding_a , forwarding_b , id_ex_rs , id_ex_rt, ex_mem_rd , mem_wb_rd , ex_mem_regwrite , mem_wb_regwrite ); 
output reg [1:0] forwarding_a; /* it will be selector in mux that selects between (reg data and last alu result and before last alu result  )*/
output reg [1:0] forwarding_b; /* it will be selector in mux that selects between (reg data and last alu result and before last alu result  )*/

input  wire [4:0] id_ex_rs; 
input  wire [4:0] id_ex_rt;
input  wire [4:0] ex_mem_rd;
input  wire [4:0] mem_wb_rd;
input  wire       ex_mem_regwrite;
input  wire       mem_wb_regwrite;

always@(id_ex_rs , id_ex_rt, ex_mem_rd , mem_wb_rd , ex_mem_regwrite , mem_wb_regwrite)
begin
/*forwarding due to 1 &2 instruction for a */
	if( (ex_mem_regwrite) && (ex_mem_rd!=0) && (ex_mem_rd==id_ex_rs) )
		begin 
			forwarding_a <= 2'b10;
		end

	else if((mem_wb_regwrite) && (mem_wb_rd !=0) && (mem_wb_rd==id_ex_rs))
		begin 
			forwarding_a <= 2'b01;
		end
	else 
		begin
			forwarding_a <= 2'b00;
		end

/*forwarding due to 1 &2 instruction for b */
	if( (ex_mem_regwrite) && (ex_mem_rd!=0) && (ex_mem_rd==id_ex_rt) )
		begin 
			forwarding_b <= 2'b10;
		end

	else if( (mem_wb_regwrite) && (mem_wb_rd !=0) && (mem_wb_rd==id_ex_rt))
		begin 
			forwarding_b <= 2'b01;
		end
	else 
		begin
			forwarding_b <= 2'b00;
		end
end 
endmodule


///////////////////////////////// pipeline processor integration module ///////////////////////
module pipeline_processor(clk);

input clk;
wire [31:0] pc;
wire [31:0] after_sign_ext_id_exe;
wire [31:0] after_sign_ext;
wire [31:0] pc_id_exe;
wire [31:0] pc_if_id;
wire [31:0] pc_plus_4;
wire taken,jump,jal,jr;
wire [31:0] inst_if_id;
wire [31:0] instruction;
wire flush;
wire mem_write;
wire reg_dst;
wire branch;
wire branch_not_eq; 
wire mem_to_reg;
wire mem_read;
wire alu_src;
wire reg_write;
wire [2:0] aluop;
wire [31:0] write_data;
wire [31:0] read_data1; //output reg_file
wire [31:0] read_data2; //output reg_file
wire [10:0] ctrl_signals_id_exe;
wire [31:0] read_data1_id_exe;
wire [31:0] read_data2_id_exe;
wire [4:0] rd_id_exe;
wire [4:0] rt_id_exe;
wire [4:0] rs_id_exe;
wire [31:0] in1; //first input to alu
wire [31:0] in2; //second input to alu
wire [31:0] out; //input to exe_mem it will use in case sw
wire [4:0]  actual_destination; // it will be rt or rd
wire [3:0]  control_bits; // output alu control
wire [31:0] alu_out;
wire zero_flag;
wire [3:0] ctrl_signals_exe_mem;
wire [4:0] reg_dst_exe_mem;
wire [31:0] alu_out_exe_mem;
wire [31:0] read_data2_exe_mem;
wire [31:0] data_mem_out;
wire [1:0] ctrl_signals_mem_wb;
wire [4:0] reg_dst_mem_wb;
wire [31:0] data_mem_out_mem_wb;
wire [31:0] alu_out_mem_wb;
wire [1:0] forwarding_a; /* it will be selector in mux that selects between (reg data and last alu result and before last alu result  )*/
wire [1:0] forwarding_b; /* it will be selector in mux that selects between (reg data and last alu result and before last alu result  )*/


pc_counter M1( pc , pc_id_exe , pc_if_id , after_sign_ext_id_exe , taken,jump , jal , jr , inst_if_id , read_data1,flush , clk );

inst_mem   M2( pc , instruction ); 

if_id_reg  M3( instruction , pc , inst_if_id , pc_if_id , taken , flush , jump , jal , jr , clk ); 

ctrl_unit  M4( inst_if_id[31:26] , mem_write , reg_dst , jump , branch , mem_to_reg , mem_read , alu_src , reg_write , aluop , branch_not_eq , jal , jr ,inst_if_id[5:0]); 

adder      M25(32'b100 , pc_if_id , pc_plus_4 );

reg_file   M5( inst_if_id[25:21] , inst_if_id[20:16] , reg_dst_mem_wb , write_data , ctrl_signals_mem_wb[0], read_data1 ,read_data2 ,jal,pc_plus_4, clk );

sign_ext   M6( inst_if_id[15:0] , after_sign_ext );


id_ex_reg  M7({ mem_write , reg_dst , branch , branch_not_eq , mem_to_reg , mem_read , alu_src , reg_write , aluop }
                , pc_if_id , read_data1 , read_data2 , after_sign_ext , inst_if_id[15:11] , inst_if_id[20:16] , inst_if_id[25:21] ,
		ctrl_signals_id_exe , pc_id_exe , read_data1_id_exe , read_data2_id_exe , after_sign_ext_id_exe , rd_id_exe ,rt_id_exe , rs_id_exe ,
		flush,taken,clk);

mux3to1    M8( read_data1_id_exe , write_data , alu_out_exe_mem , forwarding_a , in1/*first input to alu*/);

mux3to1    M9( read_data2_id_exe , write_data , alu_out_exe_mem , forwarding_b , out/*input to exe_mem it will use in case sw*/);

mux2to1   M10( out, after_sign_ext_id_exe , ctrl_signals_id_exe[4] , in2/*second input to alu*/);

mux5_2to1 M11( rt_id_exe , rd_id_exe , ctrl_signals_id_exe[9] , actual_destination );

alu_control M12( control_bits , ctrl_signals_id_exe[2:0] , after_sign_ext_id_exe [5:0] );

real_alu    M13( in1 , in2 , after_sign_ext_id_exe[10:6] , control_bits , alu_out , zero_flag );

assign taken=((zero_flag&ctrl_signals_id_exe[8])|((~(zero_flag))&ctrl_signals_id_exe[7]));

ex_mem_reg  M18( {ctrl_signals_id_exe[10], ctrl_signals_id_exe[6] , ctrl_signals_id_exe[5] , ctrl_signals_id_exe[3] }
		, alu_out , out , actual_destination
		, ctrl_signals_exe_mem , alu_out_exe_mem , read_data2_exe_mem , reg_dst_exe_mem , clk );

data_mem    M19( alu_out_exe_mem , read_data2_exe_mem , ctrl_signals_exe_mem[1] , ctrl_signals_exe_mem[3] , data_mem_out , clk);


mem_wb_reg  M20( {ctrl_signals_exe_mem[2] , ctrl_signals_exe_mem[0] }
		, reg_dst_exe_mem , data_mem_out , alu_out_exe_mem , ctrl_signals_mem_wb , reg_dst_mem_wb , data_mem_out_mem_wb , alu_out_mem_wb , clk);


mux2to1     M21( data_mem_out_mem_wb , alu_out_mem_wb , ~ctrl_signals_mem_wb[1] , write_data );


forwarding_unit M22( forwarding_a , forwarding_b , rs_id_exe , rt_id_exe , reg_dst_exe_mem , reg_dst_mem_wb , ctrl_signals_exe_mem[0] , ctrl_signals_mem_wb[0] ); 


HDU  M23(ctrl_signals_id_exe[5] , rt_id_exe , inst_if_id[25:21] , inst_if_id[20:16] , flush );

endmodule 

////////////////////// testbench just to generate the clock ////////////////////////
module pro_pipe_test;
reg clk;
pipeline_processor pp(clk);
initial
begin
$monitor("%b",clk);
clk=0;
end
always #2 clk=~clk;

endmodule
