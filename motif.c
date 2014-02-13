#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>

#include "stat.h"
#include "motif.h"

/**
  Read a Homer formatted motif and store it in memory
*/

int readMotifFile(const char *filename, Motif *m){
   FILE *mfile;
   char *line = NULL;
   size_t len = 0;
   char *saveptr;
   char *token;
   double A, T, G, C;
   int tokenlen;
   int matrixcounter;

   mfile = fopen(filename, "r");
   while( getline(&line, &len, mfile) != -1){
      if(line[0] == '>'){
         /* Process header */
         /* New entry*/
         m->count++;
         m->name   = realloc(m->name,   sizeof(char*) * m->count);
         m->length = realloc(m->length, sizeof(int)   * m->count);
         m->value  = realloc(m->value,  sizeof(char*) * m->count);

         /* new name */
         token = strtok(line+1, "\t");
         tokenlen = strlen(token);
         m->name[m->count-1] = malloc(sizeof(char) * tokenlen+1);
         strncpy(m->name[m->count-1], token, tokenlen);
         m->name[m->count-1][tokenlen] = '\0';

         /* new matrix */
         m->length[m->count-1] = tokenlen * 4;
         m->value[m->count-1] = malloc(sizeof(double) * tokenlen * 4);
         matrixcounter = 0;
      }
      else{
         /* Process body */
         token = strtok_r(line, "\t", &saveptr);
         A = atof(token);
         token = strtok_r(NULL, "\t", &saveptr);
         C = atof(token);
         token = strtok_r(NULL, "\t", &saveptr);
         G = atof(token);
         token = strtok_r(NULL, "\n", &saveptr);
         T = atof(token);

         /* Fill matrix values */
         m->value[m->count - 1][matrixcounter++] = A;
         m->value[m->count - 1][matrixcounter++] = C;
         m->value[m->count - 1][matrixcounter++] = G;
         m->value[m->count - 1][matrixcounter++] = T;
      }
   }

   free(line);
   fclose(mfile);
   return(OK);
}

Motif readDir(char *dirname){
   DIR *dir;
   struct dirent *entry;
   Motif result;
   char fullname[300];

   result.count  = 0;
   result.name   = NULL;
   result.length = NULL;
   result.value  = NULL;

   dir = opendir(dirname);
   while( (entry = readdir(dir)) != NULL){
      if(entry->d_name[0] == '.') continue;
      strcpy(fullname, dirname);
      strcat(fullname, entry->d_name);
      readMotifFile(fullname, &result);
   }

   closedir(dir);
   return(result);
}

void reverseMotif(double *a, double *b, int len){
   int i;

   for(i = 0; i < len; i++){
      b[i] = a[len-i];
   }
}

/* Compare two motifs. Make an equal length overlapping submatrix and find the best */
void pairCompare(double *a, double *b, int alen, int blen){
   int i;
   double r, logp;
   double bestr = -1e10;
   double bestlogp;
   double *reverse;

   for(i = 0; i < alen - blen + 1; i++){
      correlation(a+i, b, blen, &r, &logp);
      if(r > bestr){
         bestr    = r;
         bestlogp = logp;
      }
   }

   /* Search in reverse order */
   reverse = malloc(sizeof(double) * blen);
   reverseMotif(b, reverse, blen);

   for(i = 0; i < alen - blen + 1; i++){
      correlation(a+i, reverse, blen, &r, &logp);
      if(r > bestr){
         bestr    = r;
         bestlogp = logp;
      }
   }

   printf("\t%f\t%f", bestr, bestlogp);

   free(reverse);
}

void motifAllAgainstAll(Motif motif){
   int i, j;

   for(i = 0; i < motif.count-1; i++){
      for(j = i+1; j < motif.count; j++){
         printf("%s vs %s", motif.name[i], motif.name[j]);
         if(motif.length[i] > motif.length[j]){
            pairCompare(motif.value[i], motif.value[j], motif.length[i], motif.length[j]);
         }
         else{
            pairCompare(motif.value[j], motif.value[i], motif.length[j], motif.length[i]);
         }
         printf("\n");
      }
   }
}
