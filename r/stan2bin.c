
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define BUFLEN (64*1024*1024)
#define PARAMLEN 64

int main(int argc, char *argv[])
{
  char *cbuf;
  FILE *stan, *rdata;

  cbuf = malloc(BUFLEN);
  if (cbuf == NULL) {
    fprintf(stderr, "Memory allocation error.\n");
    return 2;
  }

  if (argc != 2) {
    printf("Usage: stan2bin STANOUTPUT\n");
    return 1;
  }

  /*
   * open stan file
   */

  snprintf(cbuf, BUFLEN, "%s.csv", argv[1]);
  stan = fopen(cbuf, "r");
  if (stan == NULL) {
    snprintf(cbuf, BUFLEN, "Error opening STAN file '%s.csv'", argv[1]);
    perror(cbuf);
    return 1;
  }

  /*
   * open binary file
   */

  snprintf(cbuf, BUFLEN, "%s.bin", argv[1]);
  rdata = fopen(cbuf, "wb");
  if (rdata == NULL) {
    snprintf(cbuf, BUFLEN, "Error opening BIN file '%s.bin'", argv[1]);
    perror(cbuf);
    return 1;
  }

  // write 3 zeros, will overwrite with proper dimensions later
  int zeros[] = { 0, 0, 0 };
  fwrite(zeros, sizeof(int), 3, rdata);

  /*
   * parse stan output
   *
   * 1. start in 'warmup' mode
   * 2. switch to 'samples' mode when we see "# Adaptation terminated"
   */

  char *params;
  float *values;
  int nparams = 0;
  int nsamples = 0;

  while (fgets(cbuf, BUFLEN, stan) != NULL) {

    if (strncmp(cbuf, "# Adapt", 7) == 0) {
      // done with warmup
      fseek(rdata, 0, SEEK_SET);
      fwrite(&nsamples, sizeof(int), 1, rdata);
      fseek(rdata, 0L, SEEK_END);
      printf("nwarmup:  %d\n", nsamples);
      nsamples = 0;
    } if (strncmp(cbuf, "lp__", 4) == 0) {
      // header line
      int ntok = 1;
      for (int i=0; i<BUFLEN; i++) {
        if (cbuf[i] == ',') ntok++;
        if (cbuf[i] == 0) break;
      }
      nparams = ntok - 2;

      printf("nparams:  %d\n", nparams);

      params = malloc(nparams*PARAMLEN);
      if (params == NULL) {
        fprintf(stderr, "Memory allocation error.\n");
        return 2;
      }

      values = malloc(nparams*sizeof(float));
      if (values == NULL) {
        fprintf(stderr, "Memory allocation error.\n");
        return 2;
      }

      char *tok;
      int i = 0;
      for (tok = strtok(cbuf, ","), ntok = 0; tok != NULL; tok = strtok(NULL, ",\n"), ntok++) {
        // skip treedepth__ and n_divergent__
        if (! (ntok == 3 || ntok == 4)) {
          strncpy(params+i*PARAMLEN, tok, PARAMLEN-1);
          i++;
        }
      }
    } else if (cbuf[0] != '#') {
      // data line
      int i = 0;
      char *tok;
      int ntok = 0;
      for (tok = strtok(cbuf, ","); tok != NULL; tok = strtok(NULL, ","), ntok++) {
        // skip treedepth__ and n_divergent__
        if (! (ntok == 3 || ntok == 4)) {
          values[i++] = (float) atof(tok);
        }
      }
      if (i == nparams) {
        fwrite(values, sizeof(float), nparams, rdata);
        nsamples++;
      }
    }
  }

  fseek(rdata, sizeof(int), SEEK_SET);
  fwrite(&nsamples, sizeof(int), 1, rdata);
  fwrite(&nparams, sizeof(int), 1, rdata);
  fseek(rdata, 0L, SEEK_END);

  printf("nsamples: %d\n", nsamples);

  for (int i=0; i<nparams; i++) {
    fwrite(params+i*PARAMLEN, sizeof(char), strlen(params+i*PARAMLEN)+1, rdata);
    //    printf("param %d: '%s'\n", i, params+i*PARAMLEN);
  }

  fclose(stan);
  fclose(rdata);

  return 0;
}
