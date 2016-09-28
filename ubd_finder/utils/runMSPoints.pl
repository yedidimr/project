#!/usr/bin/perl

if ($#ARGV != 2 && $#ARGV != 1 && $#ARGV != 0) {
  print "Usage: runMSPoints <pdb file name> [density (default-10.0)] [probe radius (default-1.8)]\n";
  exit;
}

$pdb = $ARGV[0];
if($#ARGV > 0) {
  $density = $ARGV[1];
  if($#ARGV > 1) {
    $rad = $ARGV[2];
  }
}
use File::Basename;
my $home=dirname(__FILE__);
unlink "BEFORE", "CONTACT", "REENTRANT";

if ($#ARGV == 0 ) {
  system("$home/msdots $pdb $home/vdw.lib out");
}
if ($#ARGV == 1 ) {
  system("$home/msdots $pdb $home/vdw.lib out $density");
}
if ($#ARGV == 2 ) {
  system("$home/msdots $pdb $home/vdw.lib out $density $rad");
}

open MS, ">$pdb.ms";
$number=0;
foreach $file ("CONTACT","REENTRANT") {
  open FILE, $file;
  while (<FILE>) {
    $atom1 = substr($_,0,5);
    $atom2 = substr($_,5,5);
    $atom3 = substr($_,10,5);
    $num = substr($_,15,2);
    $px = substr($_,17,9);
    $py = substr($_,26,9);
    $pz = substr($_,35,9);
    $sarea = substr($_,44,7);
    $nx = substr($_,51,7);
    $ny = substr($_,58,7);
    $nz = substr($_,65,7);
    $curv = substr($_,72,2);
    printf MS "%5d%5d%5d", $atom1, $atom2, $atom3;
    printf MS "%8.3f%8.3f%8.3f", $px, $py, $pz;
    printf MS "%8.3f", $sarea;
    printf MS "%7.3f%7.3f%7.3f", $nx, $ny, $nz;
    printf MS "%7.3f\n", 0.5;
    $number++;
  }
  close FILE;
}
close MS;

print "Surface Points number: $number\n";
unlink "BEFORE", "CONTACT", "REENTRANT";


