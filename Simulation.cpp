#include "Simulation.h"
using namespace std;

extern void read_elf(char* filename);
extern unsigned int cadr;
extern unsigned int dadr;
extern unsigned int csize;
extern unsigned long long dsize;
extern unsigned int vadr;
extern unsigned long long gp;
extern unsigned int madr;
extern unsigned int endPC;
extern unsigned int entry;
extern unsigned long long global_pointer;
extern FILE *file;

unsigned short read_mem_2(unsigned long long address);
unsigned int read_mem_4(unsigned long long address);
unsigned long long read_mem_8(unsigned long long address);
void write_mem_2(unsigned long long address, unsigned short val);
void write_mem_4(unsigned long long address, unsigned int val);
void write_mem_8(unsigned long long address, unsigned long long val);
void disp_reg();
void disp_memory(int addr, int size, int blocks);

bool print_fetch = false, print_inst = false, print_aluinfo = false, print_pc = false, print_mem = false, print_wb = false;
bool onestep = false;
//#define DBG

//指令运行数
long long inst_num=0;

//系统调用退出指示
int exit_flag=0;

//加载代码段
//初始化PC
void load_memory()
{
    fseek(file,cadr,SEEK_SET);
    fread(&memory[vadr],1,csize,file);
    fseek(file,dadr,SEEK_SET);
    fread(&memory[gp],1,dsize,file);

    //vadr=vadr>>2;
    //csize=csize>>2;
    fclose(file);
}

//run: sim -r xxx
//one-step mode: sim -r xxx -i
//print detailed information: sim -r xx -p
int main(int argc, char* argv[])
{
    char* filename = NULL;

    int oc;                     //options
    while((oc = getopt(argc, argv, "r:ip")) != -1)
    {
        switch(oc)
        {
            case 'r':
                filename = optarg;
                break;
            case 'i':
                onestep = true;
                printf("One-Step Mode.\n");
                break;
            case 'p':
                print_fetch = print_inst = print_aluinfo = print_pc = print_mem = print_wb = true;
                printf("Showing detailed information.\n");
                break;
        }
    }
    if(filename == NULL){
        cout<<"Error: No input file!"<<endl;
        return 0;
    }
    //解析elf文件
    read_elf(filename);
    //printf("%x\n", entry);
    //加载内存
    load_memory();

    //设置入口地址
    PC=entry;

    //设置全局数据段地址寄存器
    reg[3]=global_pointer;

    reg[2]=MAX/2;//栈基址 （sp寄存器）

    simulate();

    return 0;
}


bool parse_dbg_cmd(){
    //if continue, return true
    cout<<" >"<<endl;
    char dbg_cmd[100];
    cin.getline(dbg_cmd,100);
    if(dbg_cmd[0]=='c'){
        onestep = false;
        return true;
    }
    if(dbg_cmd[0]=='n'){
        return true;
    }
    if(dbg_cmd[0]=='p'){
        disp_reg();
        return false;
    }
    //x blocks size addr
    if(dbg_cmd[0]=='x'){
        int chara, size, blocks, addr;
        if(sscanf(dbg_cmd, "%c %d %d %x", &chara, &blocks, &size, &addr)){
            disp_memory(addr, size, blocks);
        }
        else{
            cout<<"Debug command error"<<endl;
        }
        return false;
    }
    cout<<"Debug command error"<<endl;
    return false;
}

void invalid_inst(){
    printf("Invalid Instruction!");
    while(true){
        bool cont = parse_dbg_cmd();
        if(cont)
            break;
    }
    exit(0);
}

void simulate()
{
    //结束PC的设置
    //int end=(int)endPC/4-1;
    int end = endPC;
    while(PC<end) {
        //运行
        IF();
        if(onestep){
            while(true){
                bool cont = parse_dbg_cmd();
                if(cont)
                    break;
            }
        }
        IF_ID = IF_ID_old;
        ID();
        ID_EX = ID_EX_old;
        EX();
        EX_MEM = EX_MEM_old;
        MEM();
        MEM_WB = MEM_WB_old;
        WB();

        //更新中间寄存器

        if (exit_flag == 1)
            break;

        reg[0] = 0;//一直为零

#ifdef DBG
        disp_reg();
        disp_memory(0x117d0, 4, 4);
#endif
        if(print_fetch || onestep) {
            printf("\n");
        }

    }
    cout <<"simulate over!"<<endl;
    while(true){
        bool cont = parse_dbg_cmd();
        if(cont)
            break;
    }
}


//取指令
void IF()
{
    //write IF_ID_old
    if(print_fetch || onestep){
        printf("Fetching instruction at %x\n", PC);
    }
    IF_ID_old.inst=read_mem_4(PC);
    IF_ID_old.PC=PC;
    PC=PC+4;
}

//译码
void ID()
{
    //Read IF_ID
    unsigned int inst=IF_ID.inst;
    int EXTop=0;
    unsigned int EXTsrc=0;

    char RegDst;
    ALUa_Source ALUSrc_a;
    ALUb_Source ALUSrc_b;
    branch_cond Branch;
    char MemRead,MemWrite;
    char RegWrite,MemtoReg,PCtoReg;
    int immlen = 0;
    ALUOPs ALUop;
    newPC_Source newPCSrc = NObranch;

    OP=getbit(inst,0,6);
    rd=getbit(inst,7,11);
    fuc3=getbit(inst,12,14);
    rs=getbit(inst,15,19);
    rt=getbit(inst,20,24);
    fuc7=getbit(inst,25,31);
    //....
    EXTsrc = 0;
    RegDst=0;
    MemRead=0;
    MemWrite=0;
    ALUSrc_a=RS;
    ALUSrc_b=RT;
    RegWrite=0;
    MemtoReg=0;
    EXTop=0;
    PCtoReg=0;
    ALUop=NONE;
    Branch=NEVER;

    if(OP==OP_R) {
        RegDst=rd;
        Branch=NEVER;
        MemRead=0;
        ALUSrc_a=RS;
        ALUSrc_b=RT;
        RegWrite=1;
        MemtoReg=0;
        EXTop=0;
        PCtoReg=0;
        if(fuc3==F3_ADD && fuc7==F7_ADD)//add
        {
            ALUop=ADD;
            if(print_inst){
                printf("Decoding: add, %d, %d, %d\n", rd,rs,rt);
            }
        }
        else if(fuc3==F3_MUL && fuc7==F7_MUL){//mul
            ALUop=MUL;
            if(print_inst){
                printf("Decoding: mul, %d, %d, %d\n", rd,rs,rt);
            }
        }
        else if(fuc3==F3_SUB && fuc7==F7_SUB){//sub
            ALUop=SUB;
            if(print_inst){
                printf("Decoding: sub, %d, %d, %d\n", rd,rs,rt);
            }
        }
        else if(fuc3==F3_SLL && fuc7==F7_SLL){//sll
            ALUop=SLL;
            if(print_inst){
                printf("Decoding: sll, %d, %d, %d\n", rd,rs,rt);
            }
        }
        else if(fuc3==F3_MULH && fuc7==F7_MULH){//mulh
            ALUop=MULH;
            if(print_inst){
                printf("Decoding: mulh, %d, %d, %d\n", rd,rs,rt);
            }
        }
        else if(fuc3==F3_SLT && fuc7==F7_SLT){//slt
            ALUop=SUB;
            if(print_inst){
                printf("Decoding: slt, %d, %d, %d\n", rd,rs,rt);
            }
        }
        else if(fuc3==F3_XOR && fuc7==F7_XOR){//xor
            ALUop=XOR;
            if(print_inst){
                printf("Decoding: xor, %d, %d, %d\n", rd,rs,rt);
            }
        }
        else if(fuc3==F3_DIV && fuc7==F7_DIV){//div
            ALUop=DIV;
            if(print_inst){
                printf("Decoding: div, %d, %d, %d\n", rd,rs,rt);
            }
        }
        else if(fuc3==F3_SRL && fuc7==F7_SRL){//srl
            ALUop=SRL;
            if(print_inst){
                printf("Decoding: srl, %d, %d, %d\n", rd,rs,rt);
            }
        }
        else if(fuc3==F3_SRA && fuc7==F7_SRA){//sra
            ALUop=SRA;
            if(print_inst){
                printf("Decoding: sra, %d, %d, %d\n", rd,rs,rt);
            }
        }
        else if(fuc3==F3_OR && fuc7==F7_OR){//or
            ALUop=OR;
            if(print_inst){
                printf("Decoding: or, %d, %d, %d\n", rd,rs,rt);
            }
        }
        else if(fuc3==F3_REM && fuc7==F7_REM){//rem
            ALUop=DIV;
            if(print_inst){
                printf("Decoding: rem, %d, %d, %d\n", rd,rs,rt);
            }
        }
        else if(fuc3==F3_AND && fuc7==F7_AND){//and
            ALUop=AND;
            if(print_inst){
                printf("Decoding: and, %d, %d, %d\n", rd,rs,rt);
            }
        }
        else
        {
            ALUop=NONE;
            invalid_inst();
        }
    }
    else if(OP==OP_I)
    {
        EXTsrc = imm12 = getbit(inst, 20, 31);
        immlen = 12;
        RegDst=rd;
        Branch=NEVER;
        MemRead=0;
        MemWrite=0;
        ALUSrc_a=RS;
        ALUSrc_b=IMM;
        RegWrite=1;
        MemtoReg=0;
        EXTop=1;
        PCtoReg=0;
        if(fuc3==F3_ADDI){//addi
            ALUop=ADD;
            if(print_inst){
                printf("Decoding: addi, %d, %d, %d\n", rd,rs,imm12);
            }
        }
        else if(fuc3==F3_SLLI && fuc7==F7_SLLI){//slli
            ALUop=SLL;
            if(print_inst){
                printf("Decoding: slli, %d, %d, %d\n", rd,rs,imm12);
            }
        }
        else if(fuc3==F3_SLTI){//slti
            ALUop=SUB;
            if(print_inst){
                printf("Decoding: slti, %d, %d, %d\n", rd,rs,imm12);
            }
        }
        else if(fuc3==F3_XORI){//xori
            ALUop=XOR;
            if(print_inst){
                printf("Decoding: xori, %d, %d, %d\n", rd,rs,imm12);
            }
        }
        else if(fuc3==F3_SRLI && fuc7==F7_SRLI){//srli
            ALUop=SRL;
            if(print_inst){
                printf("Decoding: srli, %d, %d, %d\n", rd,rs,imm12);
            }
        }
        else if(fuc3==F3_SRAI && fuc7==F7_SRAI){//srai
            ALUop=SRA;
            if(print_inst){
                printf("Decoding: srai, %d, %d, %d\n", rd,rs,imm12);
            }
        }
        else if(fuc3==F3_ORI){//ori
            ALUop=OR;
            if(print_inst){
                printf("Decoding: ori, %d, %d, %d\n", rd,rs,imm12);
            }
        }
        else if(fuc3==F3_ANDI){//andi
            ALUop=AND;
            if(print_inst){
                printf("Decoding: andi, %d, %d, %d\n", rd,rs,imm12);
            }
        }
        else{
            ALUop=NONE;
            invalid_inst();
        }
    }
    else if(OP==OP_SW)
    {
        unsigned imm4_0 = getbit(inst,7,11);
        unsigned imm11_5 = getbit(inst,25,31);
        EXTsrc = imm12 = (imm11_5<<5)|imm4_0;
        immlen = 12;
        RegDst=0;
        Branch=NEVER;
        MemRead=0;
        ALUSrc_a=RS;
        ALUSrc_b=IMM;
        RegWrite=0;
        MemtoReg=0;
        EXTop=1;
        PCtoReg=0;
        ALUop=ADD;//add
        if(fuc3==F3_SB){
            MemWrite=1;
            if(print_inst){
                printf("Decoding: sb, %d, %d(Reg %d)\n", rt, imm12, rs);
            }
        }
        else if(fuc3==F3_SH){
            MemWrite=2;
            if(print_inst){
                printf("Decoding: sh, %d, %d(Reg %d)\n", rt, imm12, rs);
            }
        }
        else if(fuc3==F3_SW){
            MemWrite=4;
            if(print_inst){
                printf("Decoding: sw, %d, %d(Reg %d)\n", rt, imm12, rs);
            }
        }
        else if(fuc3==F3_SD){
            MemWrite=8;
            if(print_inst){
                printf("Decoding: sd, %d, %d(Reg %d)\n", rt, imm12, rs);
            }
        }
        else{
            MemWrite=0;
            invalid_inst();
        }
    }
    else if(OP==OP_LW)
    {
        EXTsrc = getbit(inst, 20, 31);
        immlen = 12;
        RegDst=rd;
        Branch=NEVER;
        MemWrite=0;
        ALUSrc_a=RS;
        ALUSrc_b=IMM;
        RegWrite=1;
        MemtoReg=1;
        EXTop=1;
        PCtoReg=0;
        ALUop=ADD;//add
        if(fuc3==F3_LB){
            EXTop=1;
            MemRead=1;
            if(print_inst){
                printf("Decoding: lb, %d, %d(Reg %d)\n", rd, imm12, rs);
            }
        }
        else if(fuc3==F3_LH){
            EXTop=1;
            MemRead=2;
            if(print_inst){
                printf("Decoding: lh, %d, %d(Reg %d)\n", rd, imm12, rs);
            }
        }
        else if(fuc3==F3_LW){
            EXTop=1;
            MemRead=4;
            if(print_inst){
                printf("Decoding: lw, %d, %d(Reg %d)\n", rd, imm12, rs);
            }
        }
        else if(fuc3==F3_LD){
            MemRead=8;
            if(print_inst){
                printf("Decoding: ld, %d, %d(Reg %d)\n", rd, imm12, rs);
            }
        }
        else{
            MemRead=0;
            invalid_inst();
        }
    }
    else if(OP==OP_BEQ)
    {
        unsigned imm11_11 = getbit(inst,31,31);
        unsigned imm10_10 = getbit(inst,7,7);
        unsigned imm9_4 = getbit(inst,25,30);
        unsigned imm3_0 = getbit(inst,8,11);
        imm12 = (imm11_11<<11) | (imm10_10<<10) | (imm9_4<<4) | imm3_0;
        EXTsrc = imm12<<1;
        immlen = 13;
        RegDst=0;
        MemRead=0;
        MemWrite=0;
        ALUSrc_a=RS;
        ALUSrc_b=RT;
        RegWrite=0;
        MemtoReg=0;
        EXTop=1;
        PCtoReg=0;
        ALUop=SUB;//sub
        newPCSrc = ALUOUT;
        if(fuc3==F3_BEQ){
            Branch=EQ;
            if(print_inst){
                printf("Decoding: beq, %d, %d, %#x\n", rs, rt, imm12);
            }
        }
        else if(fuc3==F3_BNE){
            Branch=NE;
            if(print_inst){
                printf("Decoding: bne, %d, %d, %#x\n", rs, rt, imm12);
            }
        }
        else if(fuc3==F3_BLT){
            Branch=LT;
            if(print_inst){
                printf("Decoding: blt, %d, %d, %#x\n", rs, rt, imm12);
            }
        }
        else if(fuc3==F3_BGE){
            Branch=GE;
            if(print_inst){
                printf("Decoding: bge, %d, %d, %#x\n", rs, rt, imm12);
            }
        }
        else{
            invalid_inst();
            Branch=NEVER;
        }
    }
    else if(OP==OP_JAL)
    {
        unsigned imm19_19 = getbit(inst, 31, 31);
        unsigned imm9_0 = getbit(inst, 21, 30);
        unsigned imm10_10 = getbit(inst, 20, 20);
        unsigned imm18_11 = getbit(inst, 12, 19);
        imm20 = (imm19_19<<19) | (imm18_11<<11) | (imm10_10<<10) |(imm9_0);
        EXTsrc = imm20<<1;
        immlen = 21;
        RegDst=rd;
        MemRead=0;
        MemWrite=0;
        ALUSrc_a=PCval;
        ALUSrc_b=FOUR;
        RegWrite=1;
        MemtoReg=0;
        EXTop=1;
        PCtoReg=0;
        ALUop=ADD;//???
        Branch=ALWAYS;
        newPCSrc = ALUOUT;
        if(print_inst){
            printf("Decoding: jal, %d, %#x\n", rd, rt, imm20);
        }
    }
    else if(OP==OP_JALR && fuc3==F3_JALR)
    {
        imm12 = getbit(inst, 20, 31);
        EXTsrc = imm12;
        immlen = 12;
        RegDst=rd;
        MemRead=0;
        MemWrite=0;
        ALUSrc_a=PCval;
        ALUSrc_b=FOUR;
        RegWrite=1;
        MemtoReg=0;
        EXTop=1;
        PCtoReg=0;
        ALUop=ADD;//???
        Branch=ALWAYS;
        newPCSrc = ALUOUT_0;
        if(print_inst){
            printf("Decoding: jalr, %d, %d, %#x\n", rd, rs, imm12);
        }
    }
    else if(OP==OP_IW && fuc3==F3_ADDIW){
        EXTsrc = imm12 = getbit(inst, 20, 31);
        immlen = 12;
        RegDst=rd;
        Branch=NEVER;
        MemRead=0;
        MemWrite=0;
        ALUSrc_a=RS;
        ALUSrc_b=IMM;
        RegWrite=1;
        MemtoReg=0;
        EXTop=1;
        PCtoReg=0;
        ALUop=ADD;
        if(print_inst){
            printf("Decoding: addiw, %d, %d, %d\n", rd,rs,imm12);
        }
    }
    else if(OP==OP_AUIPC){
        imm20 = getbit(inst, 12, 31);
        EXTsrc = imm20<<12;
        immlen = 20;
        RegDst=rd;
        Branch=NEVER;
        MemRead=0;
        MemWrite=0;
        ALUSrc_a=PCval;
        ALUSrc_b=IMM;
        RegWrite=1;
        MemtoReg=0;
        EXTop=0;
        PCtoReg=0;
        ALUop=ADD;
        if(print_inst){
            printf("Decoding: auipc, %d, %d\n", rd,imm12);
        }
    }
    else if(OP==OP_LUI){
        imm20 = getbit(inst, 12, 31);
        EXTsrc = imm20<<12;
        immlen = 20;
        RegDst=rd;
        Branch=NEVER;
        MemRead=0;
        MemWrite=0;
        ALUSrc_a=ZERO;
        ALUSrc_b=IMM;
        RegWrite=1;
        MemtoReg=0;
        EXTop=0;
        PCtoReg=0;
        ALUop=ADD;
        if(print_inst){
            printf("Decoding: lui, %d, %#x\n", rd,imm20);
        }
    }
    else if(OP==OP_SCALL && fuc3==F3_SCALL && fuc7==F7_SCALL){
        if(reg[10]==1){
            printf("%lld\n", reg[11]);
        }
    }
    else if(OP==OP_RW){
        RegDst=rd;
        Branch=NEVER;
        MemRead=0;
        ALUSrc_a=RS;
        ALUSrc_b=RT;
        RegWrite=1;
        MemtoReg=0;
        EXTop=0;
        PCtoReg=0;
        if(fuc3==F3_ADDW && fuc7==F7_ADDW)//add
        {
            ALUop=ADD;
            if(print_inst){
                printf("Decoding: addw, %d, %d, %d\n", rd,rs,rt);
            }
        }
        else if(fuc3==F3_MULW && fuc7==F7_MULW){//mul
            ALUop=MUL;
            if(print_inst){
                printf("Decoding: mulw, %d, %d, %d\n", rd,rs,rt);
            }
        }
    }
    else{
        invalid_inst();
    }

    //write ID_EX_old
    ID_EX_old.Rd=rd;
    ID_EX_old.Rt=rt;
    ID_EX_old.Imm=ext_signed(EXTsrc,EXTop,immlen);
    ID_EX_old.Reg_Rs=reg[rs];
    ID_EX_old.Reg_Rt=reg[rt];
    ID_EX_old.PC=IF_ID.PC;
    //...

    ID_EX_old.Ctrl_EX_ALUOp=ALUop;
    ID_EX_old.Ctrl_EX_ALUSrc_a=ALUSrc_a;
    ID_EX_old.Ctrl_EX_ALUSrc_b=ALUSrc_b;
    ID_EX_old.Ctrl_EX_RegDst=RegDst;
    ID_EX_old.Ctrl_M_Branch=Branch;
    ID_EX_old.Ctrl_M_MemWrite=MemWrite;
    ID_EX_old.Ctrl_M_MemRead=MemRead;
    ID_EX_old.Ctrl_M_newPCSrc=newPCSrc;
    ID_EX_old.Ctrl_WB_MemtoReg=MemtoReg;
    ID_EX_old.Ctrl_WB_RegWrite=RegWrite;
    ID_EX_old.Ctrl_WB_PCtoReg=PCtoReg;
    //....

}

//执行
void EX()
{
    //read ID_EX
    int temp_PC=ID_EX.PC;
    char RegDst=ID_EX.Ctrl_EX_RegDst;
    char ALUOp=ID_EX.Ctrl_EX_ALUOp;
    char ALUSrc_a=ID_EX.Ctrl_EX_ALUSrc_a;
    char ALUSrc_b=ID_EX.Ctrl_EX_ALUSrc_b;

    //Branch PC calulate
    int branchPC = temp_PC + ID_EX.Imm;

    //ALU operating numbers

    REG ALU_a = 0, ALU_b = 0;

    //choose ALU input number
    //...
    if(ALUSrc_a == RS){
        ALU_a = ID_EX.Reg_Rs;
    }
    else if(ALUSrc_a == PCval){
        ALU_a = ID_EX.PC;
    }
    else if(ALUSrc_a == ZERO){
        ALU_a = 0;
    }

    if(ALUSrc_b == RT){ //using register
        ALU_b = reg[ID_EX.Rt];
    }
    else if(ALUSrc_b == IMM){//using imm
        ALU_b = (REG)ID_EX.Imm;
    }
    else if(ALUSrc_b == FOUR){
        ALU_b = 4;
    }

    //alu calculate
    int Zero = 0, Sign = 0;
    REG ALUout = 0, rem = 0;
    switch(ALUOp){
        case ADD:
            ALUout = ALU_a + ALU_b;
            break;
        case SUB:
            ALUout = ALU_a - ALU_b;
            Zero = (ALUout == 0);
            Sign = ((long long)ALU_a < (long long)ALU_b);
            break;
        case MUL:
            ALUout = ALU_a * ALU_b;
            break;
        case MULH:
            ALUout = (REG)(((__uint128_t)ALU_a*(__uint128_t)ALU_b)>>64);
        case SLL:
            ALUout = ALU_a << ALU_b;
            break;
        case SRL:
            ALUout = ALU_a >> ALU_b;
            break;
        case SRA:
            ALUout = (REG)(((long long)ALU_a)>>ALU_b);
            break;
        case AND:
            ALUout = ALU_a & ALU_b;
            break;
        case XOR:
            ALUout = ALU_a ^ ALU_b;
            break;
        case OR:
            ALUout = ALU_a | ALU_b;
            break;
        case DIV:
            ALUout = ALU_a / ALU_b;
            rem = ALU_a % ALU_b;
            break;
        default:;
    }
    if(print_aluinfo){
        printf("ALU_a = %lld, ALU_b = %lld, ALUout = %lld\n", ALU_a, ALU_b, ALUout);
    }
    //choose reg dst address
    int Reg_Dst = RegDst;

    //write EX_MEM_old
    EX_MEM_old.ALU_out=ALUout;
    EX_MEM_old.rem=rem;
    EX_MEM_old.PC=temp_PC;
    EX_MEM_old.Reg_dst=Reg_Dst;
    EX_MEM_old.Zero=Zero;
    EX_MEM_old.Sign=Sign;
    EX_MEM_old.Reg_Rt=ID_EX.Reg_Rt; //???
    EX_MEM_old.Reg_Rs=ID_EX.Reg_Rs;
    EX_MEM_old.Imm=ID_EX.Imm;
    EX_MEM_old.Ctrl_M_Branch=ID_EX.Ctrl_M_Branch;
    EX_MEM_old.Ctrl_M_MemRead=ID_EX.Ctrl_M_MemRead;
    EX_MEM_old.Ctrl_M_MemWrite=ID_EX.Ctrl_M_MemWrite;
    EX_MEM_old.Ctrl_M_newPCSrc=ID_EX.Ctrl_M_newPCSrc;
    EX_MEM_old.Ctrl_WB_RegWrite=ID_EX.Ctrl_WB_RegWrite;
    EX_MEM_old.Ctrl_WB_MemtoReg=ID_EX.Ctrl_WB_MemtoReg;
    EX_MEM_old.Ctrl_WB_PCtoReg=ID_EX.Ctrl_WB_PCtoReg;
    //.....
}

//访问存储器
void MEM()
{
    //read EX_MEM
    branch_cond Branch = EX_MEM.Ctrl_M_Branch;
    REG ALUout = EX_MEM.ALU_out;
    REG Reg_Rt = EX_MEM.Reg_Rt;
    REG Reg_Rs = EX_MEM.Reg_Rs;
    int Imm = EX_MEM.Imm;
    int Zero = EX_MEM.Zero;
    int Sign = EX_MEM.Sign;
    char MemWrite = EX_MEM.Ctrl_M_MemWrite;
    char MemRead = EX_MEM.Ctrl_M_MemRead;
    char newPCSrc = EX_MEM.Ctrl_M_newPCSrc;
    int temp_PC = EX_MEM.PC;
    int newPC = 0;
    unsigned long long Mem_read = 0;
    //complete Branch instruction PC change
    if(newPCSrc==ALUOUT){//JAL/BEQ
        newPC = temp_PC + Imm;
        //newPC = ALUout;
    }
    else if(newPCSrc==ALUOUT_0){//JALR
        newPC = (Reg_Rs + Imm) & (-2);
        //newPC = ALUout & (-2);
    }
    switch(Branch){
        case EQ:
            if(Zero) PC = newPC;
            break;
        case NE:
            if(!Zero) PC = newPC;
            break;
        case LT:
            if(Sign) PC = newPC;
            break;
        case GE:
            if(!Sign) PC = newPC;
            break;
        case ALWAYS:
            PC = newPC;
            break;
        case NEVER:
            break;
        default:;
    }
    if(print_pc){
        printf("old PC = %x, new PC = %x\n", temp_PC, PC);
    }
    //read / write memory
    unsigned long long address = ALUout;
    if(MemRead>0) {
        if (MemRead == 8) {
            Mem_read = read_mem_8(address);
        } else if (MemRead == 4) {
            Mem_read = (REG)ext_signed_64(read_mem_4(address), 1, 32);
        } else if (MemRead == 2) {
            Mem_read = (REG)ext_signed_64(read_mem_2(address), 1, 16);
        } else if (MemRead == 1) {
            Mem_read = (REG)ext_signed_64(memory[address], 1, 8);
        } else{
            printf("Memory Read Error!\n");
        }
        if(print_mem){
            printf("Read value %llx from address %llx, size = %d\n", Mem_read, address, MemRead);
        }
    }
    else if(MemWrite>0){
        if(MemWrite==8){
            write_mem_8(address, Reg_Rt);
        }
        else if(MemWrite==4){
            write_mem_4(address, (unsigned int)Reg_Rt);
        }
        else if(MemWrite==2){
            write_mem_2(address, (unsigned short)Reg_Rt);
        }
        else if(MemWrite==1){
            memory[address] = (unsigned char)Reg_Rt;
        }
        else {
            printf("Memory Write Error!\n");
        }
        if(print_mem){
            printf("Write value %llx into address %llx, size = %d\n", Reg_Rt, address, MemWrite);
        }
    }

    //write MEM_WB_old
    MEM_WB_old.Mem_read = Mem_read;
    MEM_WB_old.ALU_out = ALUout;
    MEM_WB_old.Reg_dst = EX_MEM.Reg_dst;
    MEM_WB_old.Ctrl_WB_MemtoReg = EX_MEM.Ctrl_WB_MemtoReg;
    MEM_WB_old.Ctrl_WB_RegWrite = EX_MEM.Ctrl_WB_RegWrite;
    MEM_WB_old.Ctrl_WB_PCtoReg = EX_MEM.Ctrl_WB_PCtoReg;
}


//写回
void WB()
{
    //read MEM_WB
    int Reg_Dst = MEM_WB.Reg_dst;
    char RegWrite = MEM_WB.Ctrl_WB_RegWrite;
    char MemtoReg = MEM_WB.Ctrl_WB_MemtoReg;
    char PCtoReg = MEM_WB.Ctrl_WB_PCtoReg;
    REG ALUout = MEM_WB.ALU_out, Writeval;
    unsigned long long Mem_read = MEM_WB.Mem_read;

    //write reg
    if(RegWrite){
        if(MemtoReg){
            Writeval = Mem_read;
        }
        else if(PCtoReg){
            Writeval = PC; //PC+4
        }
        else{
            Writeval = ALUout;
        }
        reg[Reg_Dst] = Writeval;
        if(print_wb){
            printf("Write val %lld into register %d\n", Writeval, Reg_Dst);
        }
    }

}

void disp_reg(){
    printf("Registers:\n");
    for(int i=0;i<32;++i){
        printf("%d: %#llx\t", i, reg[i]);
    }
    printf("\n");
}

void disp_memory(int addr, int size, int blocks){
    if(size == 1){
        for(int i=0;i<blocks;++i){
            printf("%#.8x: %#.2x\t", addr, memory[addr]);
            if((i+1)%16==0){
                printf("\n");
            }
            addr += 1;
        }
        if(blocks%16!=0){
            printf("\n");
        }
    }
    else if(size == 2){
        for(int i=0;i<blocks;++i){
            printf("%#.8x: %#.4x\t", addr, read_mem_2(addr));
            if((i+1)%8==0){
                printf("\n");
            }
            addr += 2;
        }
        if(blocks%8!=0){
            printf("\n");
        }
    }
    else if(size == 4){
        for(int i=0;i<blocks;++i){
            printf("%#.8x: %#.8x\t", addr, read_mem_4(addr));
            if((i+1)%8==0){
                printf("\n");
            }
            addr += 4;
        }
        if(blocks%8!=0){
            printf("\n");
        }
    }
    else if(size == 8){
        for(int i=0;i<blocks;++i){
            printf("%#.8x: %#.16x\t", addr, read_mem_4(addr));
            if((i+1)%4==0){
                printf("\n");
            }
            addr += 8;
        }
        if(blocks%4!=0){
            printf("\n");
        }
    }
    else{
        printf("Size Error! Only 1,2,4,8 are supported!\n");
    }
}

unsigned short read_mem_2(unsigned long long address){
    return memory[address] | (memory[address+1]<<8);
}
unsigned int read_mem_4(unsigned long long address){
    return memory[address] | (memory[address+1]<<8) | (memory[address+2]<<16) | (memory[address+3]<<24);
}
unsigned long long read_mem_8(unsigned long long address){
    return memory[address] | (memory[address+1]<<8) | (memory[address+2]<<16) | (memory[address+3]<<24)
           | (memory[address+4]<<32) | (memory[address+5]<<40) | (memory[address+6]<<48) | (memory[address+7]<<56);
}

void write_mem_2(unsigned long long address, unsigned short val){
    memory[address] = (unsigned char)val;
    memory[address+1] = (unsigned char)(val>>8);
}
void write_mem_4(unsigned long long address, unsigned int val){
    memory[address] = (unsigned char)val;
    memory[address+1] = (unsigned char)(val>>8);
    memory[address+2] = (unsigned char)(val>>16);
    memory[address+3] = (unsigned char)(val>>24);
}
void write_mem_8(unsigned long long address, unsigned long long val){
    memory[address] = (unsigned char)val;
    memory[address+1] = (unsigned char)(val>>8);
    memory[address+2] = (unsigned char)(val>>16);
    memory[address+3] = (unsigned char)(val>>24);
    memory[address+4] = (unsigned char)(val>>32);
    memory[address+5] = (unsigned char)(val>>40);
    memory[address+6] = (unsigned char)(val>>48);
    memory[address+7] = (unsigned char)(val>>56);
}