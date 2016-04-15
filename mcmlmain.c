/***********************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *  1992.
 *
 *	main program for Monte Carlo simulation of photon
 *	distribution in multi-layered turbid media.
 *
 ****/

/****
 *	THINKCPROFILER is defined to generate profiler calls in
 *	Think C. If 1, remember to turn on "Generate profiler
 *	calls" in the options menu.
 ****/
#define THINKCPROFILER 0

/* GNU cc does not support difftime() and CLOCKS_PER_SEC.*/
#define GNUCC 0

#if THINKCPROFILER
#include <profile.h>
#include <console.h>
#include <omp.h>

#endif

#include "mcml.h"

/*	Declare before they are used in main(). */
FILE *GetFile(char *);
short ReadNumRuns(FILE* );
void ReadParm(FILE* , InputStruct * );
void CheckParm(FILE* , InputStruct * );
void InitOutputData(InputStruct, OutStruct *);
void FreeData(InputStruct, OutStruct *);
double Rspecular(LayerStruct * );
void LaunchPhoton(double, LayerStruct *, PhotonStruct *);
void HopDropSpin(InputStruct  *,PhotonStruct *,tmpOutStruct *);
void SumScaleResult(InputStruct, OutStruct *);
void WriteResult(InputStruct, OutStruct, char *);
void collect(OutStruct *Out_Ptr, tmpOutStruct *cl_OUTstruct);

// Collect Function
void collect(OutStruct *Out_Ptr, tmpOutStruct *cl_OUTstruct)
{
    int i;
    if (cl_OUTstruct->Rd_valid)
        Out_Ptr->Rd_ra[cl_OUTstruct->Rdra.x][cl_OUTstruct->Rdra.y]+=cl_OUTstruct->Rdra.w;
    if (cl_OUTstruct->Tt_valid)
        Out_Ptr->Tt_ra[cl_OUTstruct->Ttra.x][cl_OUTstruct->Ttra.y]+=cl_OUTstruct->Ttra.w;        
    for (i=0; i<cl_OUTstruct->total_steps; i++)
    {
//        if (item == 43752 && (i == 702 || i == 703))
//            printf("i = %d, x = %d, y = %d, w = %E\n", i, cl_OUTstruct->data[i].x, cl_OUTstruct->data[i].y, cl_OUTstruct->data[i].w);
        Out_Ptr->A_rz[cl_OUTstruct->data[i].x][cl_OUTstruct->data[i].y]+=cl_OUTstruct->data[i].w;
    }
}

/***********************************************************
 *	If F = 0, reset the clock and return 0.
 *
 *	If F = 1, pass the user time to Msg and print Msg on
 *	screen, return the real time since F=0.
 *
 *	If F = 2, same as F=1 except no printing.
 *
 *	Note that clock() and time() return user time and real
 *	time respectively.
 *	User time is whatever the system allocates to the
 *	running of the program;
 *	real time is wall-clock time.  In a time-shared system,
 *	they need not be the same.
 *
 *	clock() only hold 16 bit integer, which is about 32768
 *	clock ticks.
 ****/
time_t PunchTime(char F, char *Msg)
{
#if GNUCC
    return(0);
#else
    static clock_t ut0;	/* user time reference. */
    static time_t  rt0;	/* real time reference. */
    double secs;
    char s[STRLEN];

    if(F==0) {
        ut0 = clock();
        rt0 = time(NULL);
        return(0);
    } else if(F==1)  {
        secs = (clock() - ut0)/(double)CLOCKS_PER_SEC;
        if (secs<0) secs=0;	/* clock() can overflow. */
        sprintf(s, "User time: %8.0lf sec = %8.2lf hr.  %s\n",
                secs, secs/3600.0, Msg);
        puts(s);
        strcpy(Msg, s);
        return(difftime(time(NULL), rt0));
    } else if(F==2) return(difftime(time(NULL), rt0));
    else return(0);
#endif
}

/***********************************************************
 *	Print the current time and the estimated finishing time.
 *
 *	P1 is the number of computed photon packets.
 *	Pt is the total number of photon packets.
 ****/
void PredictDoneTime(long P1, long Pt)
{
    time_t now, done_time;
    struct tm *date;
    char s[80];

    now = time(NULL);
    date = localtime(&now);
    strftime(s, 80, "%H:%M %x", date);
    printf("Now %s, ", s);

    done_time = now +
                (time_t) (PunchTime(2,"")/(double)P1*(Pt-P1));
    date = localtime(&done_time);
    strftime(s, 80, "%H:%M %x", date);
    printf("End %s\n", s);
}

/***********************************************************
 *	Report time and write results.
 ****/
void ReportResult(InputStruct In_Parm, OutStruct Out_Parm)
{
    char time_report[STRLEN];

    strcpy(time_report, " Simulation time of this run.");
    PunchTime(1, time_report);

    SumScaleResult(In_Parm, &Out_Parm);
    WriteResult(In_Parm, Out_Parm, time_report);
}

/***********************************************************
 *	Get the file name of the input data file from the
 *	argument to the command line.
 ****/
void GetFnameFromArgv(int argc,
                      char * argv[],
                      char * input_filename)
{
    if(argc>=2) {			/* filename in command line */
        strcpy(input_filename, argv[1]);
    } else
        input_filename[0] = '\0';
}


/***********************************************************
 *	Execute Monte Carlo simulation for one independent run.
 ****/
void DoOneRun(short NumRuns, InputStruct *In_Ptr)
{
    int i;
    long i_photon;
    /* index to photon. register for speed.*/
    OutStruct out_parm;		/* distribution of photons.*/
    PhotonStruct photon;
    long num_photons = In_Ptr->num_photons, photon_rep=10;
    tmpOutStruct* tmpOut_Ptr;
    tmpOut_Ptr = (tmpOutStruct *)malloc(sizeof(tmpOutStruct) * num_photons);
    memset(tmpOut_Ptr, 0, sizeof(tmpOutStruct) * num_photons);
    printf("4");
#if THINKCPROFILER
    InitProfile(200,200);
    cecho2file("prof.rpt",0, stdout);
#endif

    InitOutputData(*In_Ptr, &out_parm);
    //Rspecular only one time
    out_parm.Rsp = Rspecular(In_Ptr->layerspecs);
    i_photon = num_photons;
    PunchTime(0, "");

//    do {
#pragma omp parallel for private(photon)
    for(i=0 ; i<i_photon; i++)
    {
        
        /*if(num_photons - i_photon == photon_rep) {
            printf("%ld photons & %hd runs left, ", i_photon, NumRuns);
            PredictDoneTime(num_photons - i_photon, num_photons);
            photon_rep *= 10;
        }*/
        LaunchPhoton(out_parm.Rsp, In_Ptr->layerspecs, &photon);
        do    HopDropSpin(In_Ptr, &photon, &tmpOut_Ptr[i]);
        while (!photon.dead);
//    } while(--i_photon);
    }

#if THINKCPROFILER
    exit(0);
#endif
//    for (i=0; i<num_photons; i++)
//        collect(&out_parm, &tmpOut_Ptr[i]);
    ReportResult(*In_Ptr, out_parm);
    FreeData(*In_Ptr, &out_parm);
    free(tmpOut_Ptr);
}

/***********************************************************
 *	The argument to the command line is filename, if any.
 *	Macintosh does not support command line.
 ****/
char main(int argc, char *argv[])
{
    char input_filename[STRLEN];
    FILE *input_file_ptr;
    short num_runs;	/* number of independent runs. */
    InputStruct in_parm;
    int i;

    ShowVersion("Version 1.2, 1993");
    GetFnameFromArgv(argc, argv, input_filename);
    input_file_ptr = GetFile(input_filename);
    CheckParm(input_file_ptr, &in_parm);
    num_runs = ReadNumRuns(input_file_ptr);
    #pragma omp parallel 
    #pragma omp master
    {
        printf("%d threads start... \n", omp_get_num_threads());
    }
//    while(num_runs--)  {
    for( i=0; i<num_runs--; i++)
    {
        ReadParm(input_file_ptr, &in_parm);
        DoOneRun(num_runs, &in_parm);
    }

    fclose(input_file_ptr);
    return(0);
}
