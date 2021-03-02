//
//  main.c
//  genomedepth
//
//  Created by Diego Mallo on 3/26/20.
//  Copyright © 2020 DM. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <math.h>

// Expands to an integer constant expression evaluating to a close upper bound
// on the number the number of decimal digits in a value expressible in the
// integer type given by the argument (if it is a type name) or the the integer
// type of the argument (if it is an expression). The meaning of the resulting
// expression is unspecified for other arguments.
// By John Bollinger @ https://stackoverflow.com/questions/43787672/the-max-number-of-digits-in-an-int-based-on-number-of-bits
#define DECIMAL_DIGITS_BOUND(t) (241 * sizeof(t) / 100 + 1)

const int MEMBLOCK=1000000; //Size of the memory block added to the growing arrays

//My dynamic array implementation
struct di_array {
    long i;
    long allocated;
    int *mem;
};

void addvalue (struct di_array *array, int value)
{
    if (array->i>=array->allocated)
    {
        array->mem=realloc(array->mem, sizeof(*array->mem)*(array->allocated+MEMBLOCK));
        array->allocated+=MEMBLOCK;
        if (array->mem==NULL){
            EXIT_FAILURE;
        }
    }
    array->mem[array->i]=value;
    ++array->i;
}

void free_di_array (struct di_array *array){
    free(array->mem);
}
//

//My histogram implementation
struct hist_data {
    long nbins;
    double binsize;
    long min;
    long max;
    long *counts;
};

struct hist_data histogram_nbins (int *data, long ndata, long nbins)
{
    struct hist_data thehist;
    
    if (ndata==0)
    {
        thehist.counts=calloc(sizeof(*thehist.counts),nbins);
        thehist.nbins=nbins;
        thehist.max=0;
        thehist.min=0;
        thehist.binsize=0;
        return thehist;
    }
    
    int min=*data;
    int max=*(data+ndata-1);
    double bin_size=(max-min)/(double)nbins;
    thehist.counts=calloc(sizeof(*thehist.counts),nbins);
    thehist.nbins=nbins;
    thehist.max=max;
    thehist.min=min;
    thehist.binsize=bin_size;
    
    double thismax=(double)min;
    long idata=0;
    
    //[min,max) intervals
    long ibin=0;
    for (ibin=0; ibin<nbins-1; ++ibin)
    {
        thismax+=bin_size;
        while (*(data+idata)<thismax && idata<ndata)
        {
            ++*(thehist.counts+ibin);
            ++idata;
        }
    }
    
    //Last interval, [min,max]
    thismax+=bin_size;
    while (*(data+idata)<=thismax && idata<ndata)
    {
        ++*(thehist.counts+nbins-1);
        ++idata;
    }
    
    return thehist;
}

void write_histogram (struct hist_data histogram, FILE * file)
{
    long i=0;
    for (i=0; i<histogram.nbins-1;++i)
    {
        fprintf(file,"[%.2f,%.2f),%ld\n",histogram.min+histogram.binsize*i,histogram.min+histogram.binsize*(i+1),*(histogram.counts+i));
    }
    fprintf(file,"[%.2f,%.2f],%ld\n",histogram.min+histogram.binsize*(histogram.nbins-1),(double)histogram.max,*(histogram.counts+histogram.nbins-1));
}

// Compare function by balki; https://stackoverflow.com/questions/49834742/how-to-write-qsort-comparison-function-for-integers
int compare (const void *a, const void *b)
{
    const int *ia = (const int *) a;
    const int *ib = (const int *) b;
    
    if(*ia > *ib)
    {
        return 1;
    }
    else if ( *ia == *ib)
    {
        return 0;
    }
    return -1;
    
}

//My mean

double imean (int *values,long nvalues)
{
    long long total=0;
    long i=0;
    for (i=0; i<nvalues; ++i)
    {
        total+=values[i];
    }
    return total*1.0/nvalues;
}

//My total

long long itot (int *values,long nvalues)
{
    long long total=0;
    long i=0;
    for (i=0; i<nvalues; ++i)
    {
        total+=values[i];
    }
    return total;
}

//My median with pre-sorted values
double imedian (int *values, long nvalues)
{
    if (nvalues==0)
    {
        return nan("");
    }
    if (nvalues%2==0)
    {
        return (values[nvalues/2] + values [nvalues/2+1])/2.0;
    }
    else
    {
        return (double) values[nvalues/2];
    }
}

//Main program
int main(int argc, const char * argv[]) {

    char *usage="genomedepth [-i input] [-o output] [-r number_of_histogram_ranges] [-g gap_histogram_outfile] [-d depth_histogram_outfile] [-b breadth_at_depth] [-b breadth_at_depth]\nThis program parses the bedgraph file resulting from bedtools genomecov -bga. It uses stdin and stdout for input and output by default\n";
    FILE *outmessage=stderr;
    
    //Defaults
    FILE *ifile=stdin;
    FILE *ofile=stdout;
    FILE *ghfile=NULL;
    FILE *dhfile=NULL;

    int nbins=20;
    
    //ARGV parsing
    char buff_char=0; //Will contain the first character of the argument in question
    char code=0; //Argument code (i.e., i for input, o for output)
    char *ifilename=NULL;
    char *ofilename=NULL;
    char *ghfilename=NULL;
    char *dhfilename=NULL;
    int arg=0;
    int nbreadths=0;
    struct di_array *breadths=calloc(sizeof(*breadths),1);
    int parseBreadth;
    
    for (arg=1; arg<argc;++arg)
    {
        //First character of the argument
        buff_char=*argv[arg]; //Reads a character
        if (buff_char!='-' && code==0) //It is not a '-' or a known argument option
        {
            fprintf(stderr, "ERROR parsing argument %s\n%s",argv[arg],usage);
            return -1;
        }
        //New parameter
        else if (buff_char == '-') //If the argument starts with and - reads the character to recognice the option.
        {
            code=*(argv[arg]+1);
            ++arg;
            switch (toupper(code)) //If it is a lowercase converts it in the uppercase.
            {
                case 'I':
                {
                    ifilename=calloc(strlen(argv[arg])+1,sizeof(*ifilename));
                    sscanf(argv[arg],"%s",ifilename);
                    break;
                }
                case 'O':
                {
                    ofilename=calloc(strlen(argv[arg])+1,sizeof(*ofilename));
                    sscanf(argv[arg],"%s",ofilename);
                    break;
                }
                case 'D':
                {
                    dhfilename=calloc(strlen(argv[arg])+1,sizeof(*dhfilename));
                    sscanf(argv[arg],"%s",dhfilename);
                    break;
                }
                case 'G':
                {
                    ghfilename=calloc(strlen(argv[arg])+1,sizeof(*ghfilename));
                    sscanf(argv[arg],"%s",ghfilename);
                    break;
                }
                case 'R':
                {
                    sscanf(argv[arg],"%d",&nbins);
                    break;
                }
                case 'B':
                {
                    sscanf(argv[arg],"%d",&parseBreadth);
                    addvalue(breadths,parseBreadth);
                    ++nbreadths;
                    break;
                }
                default:
                {
                    fprintf(stderr, "ERROR recognizing argument %s\n",argv[arg]);
                }
                case 'H':
                {
                    fprintf(stderr, "%s",usage);
                    return -1;
                    break;
                }
            }
            code=0;
        }
    }
    
    //If the user provided input and output filenames, use them
    if (ifilename!=NULL)
    {
        ifile=fopen(ifilename, "r");
        
        if (ifile==NULL)
        {
            fprintf(stderr,"ERROR opening the input file %s\n%s",ifilename,usage);
            return -1;
        }
    }
    if (ofilename!=NULL)
    {
        ofile=fopen(ofilename, "w");
        
        if (ofile==NULL)
        {
            fprintf(stderr,"ERROR opening the output file %s\n%s",ofilename,usage);
            return -1;
        }
        
        outmessage=stdout; //If we are not writting the results to stdout, we can use it to write messages as a log
    }
    if (ghfilename!=NULL)
    {
        ghfile=fopen(ghfilename, "w");
        
        if (ghfile==NULL)
        {
            fprintf(stderr,"ERROR opening the output gap histogram file %s\n%s",ghfilename,usage);
            return -1;
        }
    }
    if (dhfilename!=NULL)
    {
        dhfile=fopen(dhfilename, "w");
        
        if (dhfile==NULL)
        {
            fprintf(stderr,"ERROR opening the output depth histogram file %s\n%s",dhfilename,usage);
            return -1;
        }
    }
    
    //Information on the run before actually starting it
    fprintf(outmessage,"\n############################################################################\nGenomedepth 1.1: calculating sequencing depth and gap statistics efficiently\n############################################################################\n\n");
    
    //Histogram information
    char ghistogram_info[1000];//Not worth spending time on making this dynamic
    if(ghfilename==NULL)
    {
        sprintf(ghistogram_info,"Gap histogram: disabled");
    }
    else
    {
        sprintf(ghistogram_info,"Gap histograms file: %s\nGap histogram ranges: %d",ghfilename,nbins);
    }
    
    //Histogram information
    char dhistogram_info[1000];//Not worth spending time on making this dynamic
    if(dhfilename==NULL)
    {
        sprintf(dhistogram_info,"Depth histogram: disabled");
    }
    else
    {
        sprintf(dhistogram_info,"Depth histograms file: %s\nDepth histogram ranges: %d",dhfilename,nbins);
    }
    
    //Breadths information
    char* breadths_info=malloc(((DECIMAL_DIGITS_BOUND(int)+2)*nbreadths+12)*sizeof(char));//Each breadth can be a maximum of int length + , and  to separate them. We add 12 characters for \tBreadths: \0
    char thisbreadth[DECIMAL_DIGITS_BOUND(int)+9]; //This should be +2, for the comma and the space. However, I am going to reuse this variable later, and I need +9 there
    
    int i=0;
    if(nbreadths==0)
    {
        sprintf(breadths_info,"Breadths: disabled");
    }
    else
    {
        sprintf(breadths_info,"Breadths: %d",*breadths->mem);
        for (i=1; i<nbreadths;++i)
        {//working here
            sprintf(thisbreadth,", %d",*(breadths->mem+i));
            strcat(breadths_info, thisbreadth);
        }
    }
    fprintf(outmessage,"Input file: %s\nOutput file: %s\n%s\n%s\n%s\n\n",ifilename==NULL?"stdin":ifilename,ofilename==NULL?"stdout":ofilename,ghistogram_info,dhistogram_info,breadths_info);
    
    //Variables for input file parsing
    ssize_t nchar=0;
    char *line=NULL;
    size_t lengthlineptr=0;
    int spos=0,fpos=0;
    long depth=0;
    char *chr=malloc(sizeof(*chr)*50);
    
    //Variables to collect input data
    struct di_array *adepths=calloc(sizeof(*adepths),1);
    struct di_array *agaps=calloc(sizeof(*agaps),1);
    long npos=0,ngaps=0,emptybases=0;
    int posthissection=0,length=0,ibreadth=0;
    long *countBreadths;
    countBreadths=calloc(nbreadths, sizeof(long));
    
    //Loop reading input file
    do
    {
        if (nchar>0) //skips empty lines
        {
            sscanf(line,"%s\t%d\t%d\t%ld\n",chr,&spos,&fpos,&depth);
            if(depth>INT_MAX)
            {
                fprintf(outmessage,"Error, the depth in this position is larger than the maximum integer number. This program needs to be re-compiled with unsigned integers or longs.");
                exit(EXIT_FAILURE);
            }
            length=fpos-spos;
            for (posthissection=0;posthissection<length;++posthissection)
            {
                addvalue(adepths,(int)depth);
                ++npos;
            }
            if(depth==0)
            {
                addvalue(agaps,length);
                ++ngaps;
                emptybases+=length;
            }
            for (ibreadth=0; ibreadth<nbreadths; ++ibreadth)
            {
                if((int)depth>=*(breadths->mem+ibreadth))
                {
                    countBreadths[ibreadth]+=length;
                }
            }
        }
        nchar=getline(&line,&lengthlineptr,ifile);

    }while(nchar!=-1);
    
    //Sorting data
    qsort(adepths->mem, npos, sizeof(*adepths->mem),&compare);
    qsort(agaps->mem, ngaps, sizeof(*agaps->mem),&compare);
    
    //Calculating mean and median
    long long totalbases=itot(adepths->mem,npos);
    double meandepth=imean(adepths->mem,npos);
    double meandepthcovered=imean(adepths->mem+emptybases,npos-emptybases);
    double meangapsize=imean(agaps->mem,ngaps);
    
    double mediandepth=imedian(adepths->mem, npos);
    double mediandepthcovered=imedian(adepths->mem+emptybases,npos-emptybases);
    double mediangapsize=imedian(agaps->mem,ngaps);
    
    //Breadths output information
    char* breadths_output_header=malloc(((DECIMAL_DIGITS_BOUND(int)+9)*nbreadths+1)*sizeof(char));//Each breadth can be a maximum of int length + ,Breadthx; +1 for the \0
    char* breadths_output_data=malloc((6*nbreadths+1)*sizeof(char));//Each breadth can be a maximum of ,100.0 + one character for \0
    
    
    if(nbreadths==0)
    {
        sprintf(breadths_output_data,"");
        sprintf(breadths_output_header,"");
    }
    else
    {
        sprintf(breadths_output_header,"Breadth%dx",*breadths->mem);
        sprintf(breadths_output_data,"%.2f",*countBreadths*100.0/npos);
        for (i=1; i<nbreadths;++i)
        {//working here
            sprintf(thisbreadth,",Breadth%dx",*(breadths->mem+i));
            strcat(breadths_output_header,thisbreadth);
            sprintf(thisbreadth,",%.2f",*(countBreadths+i)*100.0/npos);
            strcat(breadths_output_data,thisbreadth);
        }
    }
    
    fprintf(ofile,"GenomeSize,TotalBases,Breadth1x,%s,MeanDepth,MeanDepthCovered,MeanGapSize,MedianDepth,MedianDepthCovered,MedianGapSize\n%ld,%lld,%.2f,%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",breadths_output_header,npos,totalbases,100-(emptybases*100.0/npos),breadths_output_data,meandepth,meandepthcovered,meangapsize,mediandepth,mediandepthcovered,mediangapsize);
    //Calculating histograms
    
    if(ghfile!=NULL)
    {
        struct hist_data gap_hist=histogram_nbins(agaps->mem,ngaps,nbins);
        //Printing histograms
        write_histogram(gap_hist,ghfile);
    }
    
    if(dhfile!=NULL)
    {
        struct hist_data depth_hist=histogram_nbins(adepths->mem,npos,nbins);
        //Printing histograms
        write_histogram(depth_hist,dhfile);
    }

    //Freeing memory, closing files and goodbye
    free(chr);
    free(breadths_info);
    free(countBreadths);
    free(breadths_output_header);
    free(breadths_output_data);
    if (ifilename!=NULL) free(ifilename);
    if (ofilename!=NULL) free(ofilename);
    fclose(ifile);
    fclose(ofile);
    if (dhfile!=NULL) fclose (dhfile);
    if (ghfile!=NULL) fclose (ghfile);
    free_di_array(agaps);
    free_di_array(adepths);
    return 0;
}
