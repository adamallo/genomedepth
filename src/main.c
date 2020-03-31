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
    int min=*data;
    int max=*(data+ndata-1);
    double bin_size=(max-min)/(double)nbins;
    struct hist_data thehist;
    thehist.counts=calloc(sizeof(*thehist.counts),nbins);
    thehist.nbins=nbins;
    thehist.max=max;
    thehist.min=min;
    thehist.binsize=bin_size;
    
    double thismax=(double)min;
    long idata=0;
    
    //[min,max) intervals
    for (long ibin=0; ibin<nbins-1; ++ibin)
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
    for (long i=0; i<histogram.nbins-1;++i)
    {
        fprintf(file,"[%f,%f),%ld\n",histogram.min+histogram.binsize*i,histogram.min+histogram.binsize*(i+1),*(histogram.counts+i));
    }
    fprintf(file,"[%f,%f],%ld\n",histogram.min+histogram.binsize*(histogram.nbins-1),(double)histogram.max,*(histogram.counts+histogram.nbins-1));
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
    for (long i=0; i<nvalues; ++i)
    {
        total+=values[i];
    }
    return total*1.0/nvalues;
}

//My median with pre-sorted values
double imedian (int *values, long nvalues)
{
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

    char *usage="genomedepth [-i input] [-o output] [-r number_of_histogram_ranges] [-d histogram_outfile]\nThis program parses the tsv file resulting from bedtools genomecov -bga -ibam. It uses stdin and stdout for input and output by default\n";
    FILE *outmessage=stderr;
    
    //Defaults
    FILE *ifile=stdin;
    FILE *ofile=stdout;
    FILE *hfile=NULL;
    int nbins=20;
    
    //ARGV parsing
    char buff_char=0; //Will contain the first character of the argument in question
    char code=0; //Argument code (i.e., i for input, o for output)
    char *ifilename=NULL;
    char *ofilename=NULL;
    char *hfilename=NULL;
    for (int arg=1; arg<argc;++arg)
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
                    hfilename=calloc(strlen(argv[arg])+1,sizeof(*hfilename));
                    sscanf(argv[arg],"%s",hfilename);
                    break;
                }
                case 'R':
                {
                    sscanf(argv[arg],"%d",&nbins);
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
    if (hfilename!=NULL)
    {
        hfile=fopen(hfilename, "w");
        
        if (hfile==NULL)
        {
            fprintf(stderr,"ERROR opening the output histogram file %s\n%s",ofilename,usage);
            return -1;
        }
    }
    
    //Information on the run before actually starting it
    fprintf(outmessage,"Genomedepth 1.0: calculating sequencing depth and gap statistics efficiently\n");
    fprintf(outmessage,"\tInput file: %s\n\tOutput file: %s\n\tHistograms file: %s\n\tHistogram ranges: %d\n\n",ifilename==NULL?"stdin":ifilename,ofilename==NULL?"stdout":ofilename,hfilename==NULL?"Disabled":hfilename,nbins);
    
    //Variables for input file parsing
    ssize_t nchar=0;
    char *line=NULL;
    size_t lengthlineptr=0;
    int spos=0,fpos=0,depth=0;
    char *chr=malloc(sizeof(*chr)*50);
    
    //Variables to collect input data
    struct di_array *adepths=calloc(sizeof(*adepths),1);
    struct di_array *agaps=calloc(sizeof(*agaps),1);
    long npos=0,ngaps=0,emptybases=0;
    int posthissection=0,length=0;
    
    //Loop reading input file
    do
    {
        if (nchar>0) //skips empty lines
        {
            sscanf(line,"%s\t%d\t%d\t%d\n",chr,&spos,&fpos,&depth);
            length=fpos-spos;
            for (posthissection=0;posthissection<length;++posthissection)
            {
                addvalue(adepths,depth);
                ++npos;
            }
            if(depth==0)
            {
                addvalue(agaps,length);
                ++ngaps;
                emptybases+=length;
            }
        }
        nchar=getline(&line,&lengthlineptr,ifile);

    }while(nchar!=-1);
    
    //Sorting data
    qsort(adepths->mem, npos, sizeof(*adepths->mem),&compare);
    qsort(agaps->mem, ngaps, sizeof(*agaps->mem),&compare);
    
    //Calculating mean and median
    double meandepth=imean(adepths->mem,npos);
    double meandepthcovered=imean(adepths->mem+emptybases,npos-emptybases);
    double meangapsize=imean(agaps->mem,ngaps);
    
    double mediandepth=imedian(adepths->mem, npos);
    double mediandepthcovered=imedian(adepths->mem+emptybases,npos-emptybases);
    double mediangapsize=imedian(agaps->mem,ngaps);
    fprintf(ofile,"GenomeSize,Breadth1x,MeanDepth,MeanDepthCovered,MeanGapSize,MedianDepth,MedianDepthCovered,MedianGapSize\n%ld,%f,%f,%f,%f,%f,%f,%f\n",npos,emptybases*100.0/npos,meandepth,meandepthcovered,meangapsize,mediandepth,mediandepthcovered,mediangapsize);
    //Calculating histograms
    //struct hist_data depth_hist=histogram(adepths->mem,npos,nbins);
    struct hist_data gap_hist=histogram_nbins(agaps->mem,ngaps,nbins);
    //Printing histograms
    write_histogram(gap_hist,hfile);
    //Freeing memory, closing files and goodbye
    free(chr);
    if (ifilename!=NULL) free(ifilename);
    if (ofilename!=NULL) free(ofilename);
    fclose(ifile);
    fclose(ofile);
    if (hfile!=NULL) fclose (hfile);
    free_di_array(agaps);
    free_di_array(adepths);
    return 0;
}
