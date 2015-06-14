/*
 * Parse a STAN CSV file and save as a binary.  File format is:
 *
 * nwarmup    - int, 1: number of warmup samples
 * nsamples   - int, 1: number of samples
 * nparams    - int, 1: number of parameters
 * warmup     - float, nparams*nwarmup: warmup samples, row major
 * samples    - float, nparams*nsamples: samples, row major
 * parameters - char, nparams: parameter names (null terminated C strings)
 */

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define CBUFLEN  1024
#define LBUFLEN  (64*1024*1024)
#define PARAMLEN 64

int stan2bin(FILE *stan, FILE *rdata, char *lbuf)
{
  char  cbuf[CBUFLEN];
  int   i;
  int   ntok;
  char *tok;
  int   zeros[] = { 0, 0, 0 };

  char *params = NULL;
  float *values = NULL;
  int nparams = 0;
  int nwarmup = 0;
  int nsamples = 0;

  /*
   * write 3 zeros (counts), will overwrite with proper dimensions later
   */
  fwrite(zeros, sizeof(int), 3, rdata);

  /*
   * parse stan output
   *
   * 1. start in 'warmup' mode
   * 2. switch to 'samples' mode when we see "# Adaptation terminated"
   */

  while (fgets(lbuf, LBUFLEN, stan) != NULL) {

    if (strncmp(lbuf, "lp__", 4) == 0) {

      /*
       * header line, count parameters
       */

      for (i=0, ntok=1; i<LBUFLEN; i++) {
        if (lbuf[i] == ',') ntok++;
        if (lbuf[i] ==   0) break;
      }
      nparams = ntok;

      printf("nparams:  %d\n", nparams);

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

    } else if (strncmp(lbuf, "# Adapt", 7) == 0) {

      /*
       * end of warmup
       */

      nwarmup = nsamples;
      nsamples = 0;

      printf("nwarmup:  %d\n", nwarmup);

    } else if ((lbuf[0] != '#') && isgraph(lbuf[0])) {

      /*
       * parse sample and write
       */
      for (ntok=0, tok=strtok(lbuf, ","); tok != NULL; ntok++, tok=strtok(NULL, ",\n")) {
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
   * write parameter names and dimensions
   */

  for (i=0; i<nparams; i++) {
    fwrite(params+i*PARAMLEN, sizeof(char), strlen(params+i*PARAMLEN)+1, rdata);
  }

  for (i=0; i<nparams; i++) {
    tok = strtok(params+i*PARAMLEN, ".");
    fwrite(tok, sizeof(char), strlen(tok)+1, rdata);
  }

  fseek(rdata, 0L, SEEK_SET);
  fwrite(&nwarmup, sizeof(int), 1, rdata);
  fwrite(&nsamples, sizeof(int), 1, rdata);
  fwrite(&nparams, sizeof(int), 1, rdata);

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
