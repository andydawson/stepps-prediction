
read_stanbin <- function(fname) {
  bin = file(fname, "rb")
  nwarmup  = readBin(bin, "integer")
  nsamples = readBin(bin, "integer")
  nparams  = readBin(bin, "integer")
  warmup   = readBin(bin, "numeric", n=nwarmup*nparams, size=4)
  samples  = readBin(bin, "numeric", n=nsamples*nparams, size=4)
  params   = readBin(bin, "character", n=nparams)
  short    = readBin(bin, "character", n=nparams)
  close(bin)

  dim(warmup) <- c(nparams, nwarmup)
  dim(samples) <- c(nparams, nsamples)

  rownames(warmup) <- params
  rownames(samples) <- params

  list(params=params, warmup=t(warmup), samples=t(samples), par_names=short)
}
