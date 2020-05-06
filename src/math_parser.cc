#include "math_parser.h"

#include <cassert>
#include <iostream>
#include <string>

#include "constants.h"
#include "muParserDLL.h"
#include "ui.h"

#define PARSER_MAXVARS		500

#ifndef _UNICODE
    #define _T(x) x
    #define myprintf printf
    #define mystrlen strlen
    #define myfgets fgets
    #define mystrcmp strcmp
#else
#define _T(x) L ##x
    #define myprintf wprintf
    #define mystrlen wcslen
    #define myfgets fgetws
    #define mystrcmp wcscmp
#endif

//extern double vl_emit_nx, vl_emit_ny, vl_dp_p, vl_sigma_s, vl_rx_ibs, vl_ry_ibs, vl_rs_ibs,
//    vl_rx_ecool, vl_ry_ecool, vl_rs_ecool, vl_rx_total, vl_ry_total, vl_rs_total, vl_t;

extern Record uircd;

//---------------------------------------------------------------------------
// Factory function for creating new parser variables
// This could as well be a function performing database queries.
muFloat_t* AddVariable(const muChar_t* a_szName, void *pUserData)
{
    static muFloat_t afValBuf[PARSER_MAXVARS];  // I don't want dynamic allocation here
    static int iVal = 0;                     // so i used this buffer

//    myprintf(_T("Generating new variable \"%s\" (slots left: %d; context pointer: 0x%x)\n"), a_szName, PARSER_MAXVARS - iVal, (int)pUserData);
//    myprintf(_T("Generating new variable \"%s\" (slots left: %d \n"), a_szName, PARSER_MAXVARS - iVal);

    afValBuf[iVal] = 0;
    if (iVal >= PARSER_MAXVARS - 1)
    {
        myprintf(_T("Variable buffer overflow."));
        return NULL;
    }

    return &afValBuf[iVal++];
}

//---------------------------------------------------------------------------
void ListVar(muParserHandle_t a_hParser)
{
    int iNumVar = mupGetVarNum(a_hParser);
    int i = 0;

    if (iNumVar == 0)
    {
        myprintf(_T("No variables defined\n"));
        return;
    }

    myprintf(_T("\nExpression variables:\n"));
    myprintf(_T("---------------------\n"));
    myprintf(_T("Number: %d\n"), iNumVar);

    for (i = 0; i < iNumVar; ++i)
    {
        const muChar_t* szName = 0;
        muFloat_t* pVar = 0;

        mupGetVar(a_hParser, i, &szName, &pVar);
        myprintf(_T("Name: %s    Address: [0x%llx]\n"), szName, (unsigned long long int)pVar);
    }
}

muFloat_t Sin(muFloat_t v){return sin(v);}
muFloat_t Cos(muFloat_t v){return cos(v);}
muFloat_t Tan(muFloat_t v){return tan(v);}
muFloat_t Asin(muFloat_t v){return asin(v);}
muFloat_t Acos(muFloat_t v){return acos(v);}
muFloat_t Atan(muFloat_t v){return atan(v);}
muFloat_t Sinh(muFloat_t v){return sinh(v);}
muFloat_t Cosh(muFloat_t v){return cosh(v);}
muFloat_t Tanh(muFloat_t v){return tanh(v);}
muFloat_t Asinh(muFloat_t v){return asinh(v);}
muFloat_t Acosh(muFloat_t v){return acosh(v);}
muFloat_t Atanh(muFloat_t v){return atanh(v);}
muFloat_t Log(muFloat_t v){return log(v);}
muFloat_t Log10(muFloat_t v){return log10(v);}
muFloat_t Log2(muFloat_t v){return log2(v);}
muFloat_t Exp(muFloat_t v){return exp(v);}
muFloat_t Sqrt(muFloat_t v){return sqrt(v);}
muFloat_t Sign(muFloat_t v){if(v>0) return 1; if(v<0) return -1; return 0;}
muFloat_t Rint(muFloat_t v){return rint(v);}
muFloat_t Floor(muFloat_t v){return floor(v);}
muFloat_t Ceil(muFloat_t v){return ceil(v);}
muFloat_t Abs(muFloat_t v){return fabs(v);}

muFloat_t Sum(const muFloat_t *a_afArg, int a_iArgc){
    muFloat_t fRes = 0;
    int i = 0;
    for (i = 0; i < a_iArgc; ++i)
        fRes += a_afArg[i];
    return fRes;
}

muFloat_t Avg(const muFloat_t *a_afArg, int a_iArgc){
    muFloat_t fRes = 0;
    int i = 0;
    for (i = 0; i < a_iArgc; ++i)
        fRes += a_afArg[i];
    fRes /= a_iArgc;
    return fRes;
}

muFloat_t Max(const muFloat_t *a_afArg, int a_iArgc){
    muFloat_t fRes = a_afArg[0];
    int i = 1;
    for (i = 1; i < a_iArgc; ++i)
        if(a_afArg[i]>fRes)
            fRes = a_afArg[i];
    return fRes;
}

muFloat_t Min(const muFloat_t *a_afArg, int a_iArgc){
    muFloat_t fRes = a_afArg[0];
    int i = 1;
    for (i = 1; i < a_iArgc; ++i)
        if(a_afArg[i]<fRes)
            fRes = a_afArg[i];
    return fRes;
}

//---------------------------------------------------------------------------
void ListExprVar(muParserHandle_t a_hParser)
{
    muInt_t iNumVar = mupGetExprVarNum(a_hParser),
        i = 0;

    if (iNumVar == 0)
    {
        myprintf(_T("Expression dos not contain variables\n"));
        return;
    }

    myprintf(_T("\nExpression variables:\n"));
    myprintf(_T("---------------------\n"));
    myprintf(_T("Expression: %s\n"), mupGetExpr(a_hParser));
    myprintf(_T("Number: %d\n"), iNumVar);

    for (i = 0; i < iNumVar; ++i)
    {
        const muChar_t* szName = 0;
        muFloat_t* pVar = 0;

        mupGetExprVar(a_hParser, i, &szName, &pVar);
        myprintf(_T("Name: %s   Address: [0x%llx]\n"), szName, (unsigned long long int)pVar);
    }
}

//---------------------------------------------------------------------------
void ListConst(muParserHandle_t a_hParser)
{
    muInt_t iNumVar = mupGetConstNum(a_hParser),
        i = 0;

    if (iNumVar == 0)
    {
        myprintf(_T("No constants defined\n"));
        return;
    }

    myprintf(_T("\nParser constants:\n"));
    myprintf(_T("---------------------\n"));
    myprintf(_T("Number: %d\n"), iNumVar);

    for (i = 0; i < iNumVar; ++i)
    {
        const muChar_t* szName = 0;
        muFloat_t fVal = 0;

        mupGetConst(a_hParser, i, &szName, &fVal);
        myprintf(_T("  %s = %e\n"), szName, fVal);
    }
}



//---------------------------------------------------------------------------
// Callback function for parser errors
void OnError(muParserHandle_t hParser)
{
    myprintf(_T("\nError:\n"));
    myprintf(_T("------\n"));
    myprintf(_T("Message:  \"%s\"\n"), mupGetErrorMsg(hParser));
    myprintf(_T("Token:    \"%s\"\n"), mupGetErrorToken(hParser));
    myprintf(_T("Position: %d\n"), mupGetErrorPos(hParser));
    myprintf(_T("Errc:     %d\n"), mupGetErrorCode(hParser));
    assert(false&&"ERROR IN MATH PARSING!");
}

void AddConst(muParserHandle_t hParser) {
    mupDefineConst(hParser, "K_C", k_c);
    mupDefineConst(hParser, "K_E", k_e);
    mupDefineConst(hParser, "K_PI", k_pi);
    mupDefineConst(hParser, "K_U", k_u);
    mupDefineConst(hParser, "K_ME", k_me);
    mupDefineConst(hParser, "K_KE", k_ke);
    mupDefineConst(hParser, "k_c", k_c);
    mupDefineConst(hParser, "k_e", k_e);
    mupDefineConst(hParser, "k_pi", k_pi);
    mupDefineConst(hParser, "k_u", k_u);
    mupDefineConst(hParser, "k_me", k_me);
    mupDefineConst(hParser, "k_ke", k_ke);
}

void AddFun(muParserHandle_t hParser) {
    mupDefineFun1(hParser, _T("SIN"), Sin, 1);
    mupDefineFun1(hParser, _T("COS"), Cos, 1);
    mupDefineFun1(hParser, _T("TAN"), Tan, 1);
    mupDefineFun1(hParser, _T("ASIN"), Asin, 1);
    mupDefineFun1(hParser, _T("ACOS"), Acos, 1);
    mupDefineFun1(hParser, _T("ATAN"), Atan, 1);
    mupDefineFun1(hParser, _T("SINH"), Sinh, 1);
    mupDefineFun1(hParser, _T("COSH"), Cosh, 1);
    mupDefineFun1(hParser, _T("TANH"), Tanh, 1);
    mupDefineFun1(hParser, _T("ASINH"), Asinh, 1);
    mupDefineFun1(hParser, _T("ACOSH"), Acosh, 1);
    mupDefineFun1(hParser, _T("ATANH"), Atanh, 1);
    mupDefineFun1(hParser, _T("Ln"), Log, 1);
    mupDefineFun1(hParser, _T("LOG"), Log, 1);
    mupDefineFun1(hParser, _T("LOG10"), Log10, 1);
    mupDefineFun1(hParser, _T("LOG2"), Log2, 1);
    mupDefineFun1(hParser, _T("EXP"), Exp, 1);
    mupDefineFun1(hParser, _T("SQRT"), Sqrt, 1);
    mupDefineFun1(hParser, _T("RINT"), Rint, 1);
    mupDefineFun1(hParser, _T("FLOOR"), Floor, 1);
    mupDefineFun1(hParser, _T("CEIL"), Ceil, 1);
    mupDefineFun1(hParser, _T("floor"), Floor, 1);
    mupDefineFun1(hParser, _T("ceil"), Ceil, 1);
    mupDefineFun1(hParser, _T("ABS"), Abs, 1);
    mupDefineMultFun(hParser, _T("SUM"), Sum, 1);
    mupDefineMultFun(hParser, _T("AVG"), Avg, 1);
    mupDefineMultFun(hParser, _T("MEAN"), Avg, 1);
    mupDefineMultFun(hParser, _T("MAX"), Max, 1);
    mupDefineMultFun(hParser, _T("MIN"), Min, 1);
}


//Initialize the math parser environment
void initialize_parser(muParserHandle_t &math_parser) {
    math_parser = mupCreate(muBASETYPE_FLOAT);              // initialize the parser
    mupSetErrorHandler(math_parser, OnError);
    mupSetVarFactory(math_parser, AddVariable, NULL);       // Set a variable factory
    AddConst(math_parser);
    AddFun(math_parser);                                // Add constants

    //Set variables to save the calculation/simulation results.
    mupDefineVar(math_parser, "VL_EMIT_NX", &uircd.emit_nx);
    mupDefineVar(math_parser, "VL_EMIT_NY", &uircd.emit_ny);
    mupDefineVar(math_parser, "VL_MOMENTUM_SPREAD", &uircd.dp_p);
    mupDefineVar(math_parser, "VL_BUNCH_LENGTH", &uircd.sigma_s);
    mupDefineVar(math_parser, "VL_RATE_IBS_X", &uircd.rx_ibs);
    mupDefineVar(math_parser, "VL_RATE_IBS_Y", &uircd.ry_ibs);
    mupDefineVar(math_parser, "VL_RATE_IBS_S", &uircd.rs_ibs);
    mupDefineVar(math_parser, "VL_RATE_ECOOL_X", &uircd.rx_ecool);
    mupDefineVar(math_parser, "VL_RATE_ECOOL_Y", &uircd.ry_ecool);
    mupDefineVar(math_parser, "VL_RATE_ECOOL_S", &uircd.rs_ecool);
    mupDefineVar(math_parser, "VL_RATE_TOTAL_X", &uircd.rx_total);
    mupDefineVar(math_parser, "VL_RATE_TOTAL_Y", &uircd.ry_total);
    mupDefineVar(math_parser, "VL_RATE_TOTAL_S", &uircd.rs_total);
    mupDefineVar(math_parser, "VL_T", &uircd.t);
}

int use_parser(int argc, char* argv[])
{
    muChar_t szLine[100];
    muFloat_t fVal = 0;
    muParserHandle_t hParser;
    hParser = mupCreate(muBASETYPE_FLOAT);              // initialize the parser
    mupSetErrorHandler(hParser, OnError);

    // Set a variable factory
    mupSetVarFactory(hParser, AddVariable, NULL);

    AddConst(hParser);


    ListConst(hParser);

    std::string str = "x = 10.6e-1";
    mupSetExpr(hParser, str.c_str());

    fVal = mupEval(hParser);

    ListVar(hParser);
    ListExprVar(hParser);

    return 0;
}
