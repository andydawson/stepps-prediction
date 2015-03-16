#!/usr/bin/env perl

open(STAN, "+<", $ARGV[0]) or die;

# count samples and figure out if we need timing info
$samples = 0;
$needs_timing = 1;
while (<STAN>) {
    $counting = 1 if /Adaptation terminated/;
    $counting = 0 if /Elapsed Time/;
    $count++ if $counting;
    $needs_timing = 0 if /Elapsed Time/;
}
$samples = $count - 4;

if ($needs_timing) {
    print STAN "#  Elapsed Time: 666 seconds (Warm-up)\n";
    print STAN "#                666 seconds (Sampling)\n";
    print STAN "#                666 seconds (Total)\n";

# rewind and replace "num_samples"
    seek(STAN, 0, SEEK_SET);
    while (<STAN>) {
	$len = length();
	if (/num_samples = /) {
	    seek(STAN, tell(STAN)-$len, SEEK_SET);
	    $new = sprintf("#     num_samples = %d", $samples);
	    $new = sprintf("%-*s\n", $len-1, $new);
	    print STAN $new;
	    last;
	}
    }
}

close(STAN);
