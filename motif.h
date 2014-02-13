typedef struct _Motif{
   char   **name;     /* names of the motifs */
   int     *length;   /* length of the matrix of the motifs */
   double **value;    /* matrix of motifs */
   int      count;    /* Number of motifs */
} Motif;

Motif readDir(char *dirname);
void motifAllAgainstAll(Motif motif);
