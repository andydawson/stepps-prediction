/*
 * Parse a STAN CSV file and save as a binary.  File format is:
 *
 * magic      - int, 1: negative version counter
 * nwarmup    - int, 1: number of warmup samples
 * nsamples   - int, 1: number of samples
 * nparams    - int, 1: number of parameters
 * ndiag      - int, 1: number of diagonal entries
 * warmup     - float, nparams*nwarmup: warmup samples, row major
 * step size  - float: step size
 * diagonal   - float, ndiag: diagonal entries of inverse mass matrix 
 * samples    - float, nparams*nsamples: samples, row major
 * parameters - char, nparams: parameter names (null terminated C strings)
 */

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define CBUFLEN  1024
#define LBUFLEN  (256*1024*1024)
#define PARAMLEN 64

int stan2bin(FILE *stan, FILE *rdata, char *lbuf)
{
  char  cbuf[CBUFLEN];
  int   i;
  int   ntok;
  char *tok;
  int   header[] = { -1, 0, 0, 0, 0 };
  int   nline;

  char *params = NULL;
  float *values = NULL;
  int nparams = 0;
  int nwarmup = 0;
  int nsamples = 0;
  int ndiag = 0;
  int diag = 0;

  /*
   * write header (counts), will overwrite zeros with proper dimensions later
   */
  fwrite(header, sizeof(int), 5, rdata);

  /*
   * parse stan output
   *
   * 1. start in 'warmup' mode
   * 2. switch to 'samples' mode when we see "# Adaptation terminated"
   */

  nline = 0;
  while (fgets(lbuf, LBUFLEN, stan) != NULL) {
    nline++;

    if (strlen(lbuf) <= 0) continue;

    if (strncmp(lbuf, "lp__", 4) == 0) {

      /*
       * header line, count parameters
       */

      for (i=0, ntok=1; i<LBUFLEN; i++) {
        if (lbuf[i] == ',') ntok++;
        if (lbuf[i] ==   0) break;
      }
      nparams = ntok;

      printf("nparams:  %d (%d)\n", nparams-6, nparams);

      // allocate buffers
      params = malloc(nparams*PARAMLEN);
      if (params == NULL) {
        snprintf(cbuf, CBUFLEN,
                 "Memory allocation error (params, %d bytes)", nparams*PARAMLEN);
        perror(cbuf);
        return 2;
      }

      values = malloc(nparams*sizeof(float));
      if (values == NULL) {
        snprintf(cbuf, CBUFLEN,
                 "Memory allocation error (values, %ld bytes)", nparams*sizeof(float));
        perror(cbuf);
        return 2;
      }

      // copy parameter names
      for (ntok=0, tok=strtok(lbuf, ","); tok != NULL; ntok++, tok=strtok(NULL, ",\n")) {
	strncpy(params+ntok*PARAMLEN, tok, PARAMLEN-1);
      }

    } else if (diag == 1) {

      /* 
       * parse diagonal and save
       */
      for (ntok=0, tok=strtok(lbuf+2, ","); tok != NULL; ntok++, tok=strtok(NULL, ",\n")) {
	if (ntok >= nparams) {
	  fprintf(stderr, "Skipping excess diagonal entries on line %d of STAN output file.\n", nline);
	  break;
	}
	values[ntok] = (float) atof(tok);
      }
      ndiag = ntok;

      fwrite(values, sizeof(float), ndiag, rdata);

      diag = 0;

    } else if (strncmp(lbuf, "# Diago", 7) == 0) {

      diag = 1;
      
    } else if (strncmp(lbuf, "# Step size = ", 14) == 0) {

      float step_size = atof(lbuf+14);
      fwrite(&step_size, sizeof(float), 1, rdata);

    } else if (strncmp(lbuf, "# Adapt", 7) == 0) {

      /*
       * end of warmup
       */

      nwarmup = nsamples;
      nsamples = 0;

      printf("nwarmup:  %d\n", nwarmup);

    } else if (strncmp(lbuf, "#  Elapsed", 10) == 0) {

      break;

    } else if ((lbuf[0] != '#') && isgraph(lbuf[0])) {

      /*
       * parse sample and write
       */
      for (ntok=0, tok=strtok(lbuf, ","); tok != NULL; ntok++, tok=strtok(NULL, ",\n")) {
	if (ntok >= nparams) {
	  //	  fprintf(stderr, "ERROR: Excess parameters on line %d of output file (%d/%d).\n", nline, ntok, nparams);
	  fprintf(stderr, "Skipping excess parameters on line %d of STAN output file.\n", nline);
	  break;
	}
	values[ntok] = (float) atof(tok);
      }
      if (ntok == nparams) {
        fwrite(values, sizeof(float), nparams, rdata);
        nsamples++;
      }

    }

  }

  printf("nsamples: %d\n", nsamples);

  /*
   * write parameter names, diagonal values, and dimensions
   */

  for (i=0; i<nparams; i++) {
    fwrite(params+i*PARAMLEN, sizeof(char), strlen(params+i*PARAMLEN)+1, rdata);
  }

  for (i=0; i<nparams; i++) {
    tok = strtok(params+i*PARAMLEN, ".");
    fwrite(tok, sizeof(char), strlen(tok)+1, rdata);
  }

  fseek(rdata, sizeof(int), SEEK_SET);
  fwrite(&nwarmup, sizeof(int), 1, rdata);
  fwrite(&nsamples, sizeof(int), 1, rdata);
  fwrite(&nparams, sizeof(int), 1, rdata);
  fwrite(&ndiag, sizeof(int), 1, rdata);

  free(values);
  free(params);

  return 0;
}


int main(int argc, char *argv[])
{
  char cbuf[CBUFLEN];
  char *lbuf;
  FILE *stan, *rdata;

  if (argc != 2) {
    printf("Usage: stan2bin STANOUTPUT\n");
    return 1;
  }

  /*
   * open stan and bin files
   */

  snprintf(cbuf, CBUFLEN, "%s.csv", argv[1]);
  stan = fopen(cbuf, "r");
  if (stan == NULL) {
    snprintf(cbuf, CBUFLEN, "Error opening STAN file '%s.csv'", argv[1]);
    perror(cbuf);
    return 2;
  }

  snprintf(cbuf, CBUFLEN, "%s.bin", argv[1]);
  rdata = fopen(cbuf, "wb");
  if (rdata == NULL) {
    snprintf(cbuf, CBUFLEN, "Error opening BIN file '%s.bin'", argv[1]);
    perror(cbuf);
    return 2;
  }

  /*
   * allocate line buffer
   */

  lbuf = malloc(LBUFLEN);
  if (lbuf == NULL) {
    snprintf(cbuf, CBUFLEN, "Memory allocation error (line buffer, %d bytes)", LBUFLEN);
    perror(cbuf);
    return 3;
  }

  /*
   * parse and tidy up
   */

  int r = stan2bin(stan, rdata, lbuf);

  free(lbuf);
  fclose(stan);
  fclose(rdata);

  return r;
}
