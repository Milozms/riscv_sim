typedef unsigned long long REG;
enum ALUOPs {NONE, ADD, SUB, MUL, MULH, SLL, SRA, SRL, AND, XOR, OR, DIV};
enum branch_cond {NEVER, ALWAYS, EQ, NE, LT, GE};
enum newPC_Source {NObranch, ALUOUT, ALUOUT_0};
enum ALUa_Source {RS, PCval, ZERO};
enum ALUb_Source {RT, IMM, FOUR};

struct IFID{
	unsigned int inst;
	int PC;
}IF_ID,IF_ID_old;


struct IDEX{
	int Rd,Rt;
	int PC;
	int Imm;
	REG Reg_Rs,Reg_Rt;

	ALUa_Source Ctrl_EX_ALUSrc_a;
    ALUb_Source Ctrl_EX_ALUSrc_b;
	ALUOPs Ctrl_EX_ALUOp;
	char Ctrl_EX_RegDst;

	branch_cond Ctrl_M_Branch; //condition for branch
	char Ctrl_M_MemWrite; //num of bytns written
	char Ctrl_M_MemRead;  //num of bytes read
    newPC_Source Ctrl_M_newPCSrc;

	char Ctrl_WB_RegWrite;
	char Ctrl_WB_MemtoReg;
    char Ctrl_WB_PCtoReg;
    int target;

}ID_EX,ID_EX_old;

struct EXMEM{
	int PC;
	int Reg_dst, Imm;
	REG ALU_out, rem;
	int Zero, Sign;
	REG Reg_Rt, Reg_Rs;

	branch_cond Ctrl_M_Branch; //condition for branch
	char Ctrl_M_MemWrite; //num of bytns written
	char Ctrl_M_MemRead;  //num of bytes read
    newPC_Source Ctrl_M_newPCSrc;

	char Ctrl_WB_RegWrite;
	char Ctrl_WB_MemtoReg;
    char Ctrl_WB_PCtoReg;

}EX_MEM,EX_MEM_old;

struct MEMWB{
	unsigned long long Mem_read;
	REG ALU_out;
	int Reg_dst;
		
	char Ctrl_WB_RegWrite;
	char Ctrl_WB_MemtoReg;
    char Ctrl_WB_PCtoReg;

}MEM_WB,MEM_WB_old;