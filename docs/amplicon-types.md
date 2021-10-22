# Amplicon Types

AmpliGone requires user information when it comes to the 'amplicon type'. This describes the length of a PCR amplicon, and which part of said amplicon is actually sequenced. This information is required because it tells AmpliGone from which side(s) of a read the primer-sequence has to be removed.

The given options are "End-to-End" and "End-to-Mid".

## End-to-End

An "End-to-End" amplicon type describes the situation where a read covers the full length of an amplicon. This is most prominent in, but not exclusive to, long read sequencing platforms such as Nanopore.  
In this case of an "End-to-End" amplicon type, primers need to be removed from both ends of a read.

![](img/amplicon-type_end-to-end.png)

When given the "End-to-End" amplicon type, AmpliGone follows the logic that both forward and reverse reads should start at a forward-primer and stop at a reverse-primer. The orientation of the reads doesn't matter in this case.

## End-to-Mid

An "End-to-Mid" amplicon type describes the situation where a reads cover the amplicon only partially, and partially overlap in the amplicon mid-point. This is most prominent in, but not exclusive to, short read Illumina sequencing.  
In this case of an "End-to-Mid" amplicon type, primers need to be removed from only one end of a read.

![](img/amplicon-type_end-to-mid.png)

When given the "End-to-Mid" amplicon type, AmpliGone follows the logic that a forward-read may only start at a forward-primer. And a reverse-read may only start at at a reverse-primer.

## Why is a more common terminology not used?

Because we quickly found out that more common terminology such as paired-end/single-end, or short/long read, doesn't *always* convey the information that AmpliGone needs in a single keyword or argument. And information such as a sequencing platform doesn't always convey this information either.

As an example: it's most common that primers are only present on one side of a read if the sequencing was performed with "short reads" on an Illumina platform. ***But that's not necessarily always the case***, because depending on the actual experiment setup the primers might actually be on both ends of those short reads.

With AmpliGone we wanted a way to actually get the information that is required and matches with an experiment design, without the user having to give multiple pieces of information through various command-line arguments.