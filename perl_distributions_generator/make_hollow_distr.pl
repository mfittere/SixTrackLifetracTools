# Script for generating hollow beam distribution 
# in Lifetrac Norm format
#
use Math::Trig;
use Math::Random;
use Switch;
use GD::Graph::histogram;
use POSIX qw/strftime/;
#
if( $ARGV[0] eq "" ){
    print "usage: perl make_hollow_distr.pl distr_out\n ";
    exit(0);
}
open(fpw, ">".$ARGV[0]) || die "Cannot open $ARGV[0] $!\n";
#
$r1=4;       # initial radius (sigma)
$r2=6.7443;  # final radius
$np=10000;   # number of particles
#$zd="gauss"; # type of longitudinal distribution
             # "gauss" for Gaussian with sigma=1
$zd="zero";  # "zero" - no length and momentum spread
#$zd="delta1";# "delta1" - delta function at 1 sigma
#$zd="delta2";# "delta2" - delta function at 2 sigma
#$zd="delta3";# "delta3" - ----------------- 3 sigma
#

printf fpw "Norm\n";

random_set_seed((1,1));
@za=random_normal($np,0,1);
@zp=random_normal($np,0,1);
@rr=random_uniform($np,$r1,$r2);
@rf=random_uniform($np);

for($i=1;$i<=$np;$i++){
 switch ( $zd ){
     case "gauss" {
	 $z=$za[$i-1];
	 $pz=$zp[$i-1];
	 $x=$rr[$i-1]*cos(2*pi*$rf[$i-1]);
	 $y=$rr[$i-1]*sin(2*pi*$rf[$i-1]);
	 $px=0;
	 $py=0;
     }
     case "zero" {
	 $z=0;
	 $pz=0;
	 $x=$rr[$i-1]*cos(2*pi*$rf[$i-1]);
	 $y=$rr[$i-1]*sin(2*pi*$rf[$i-1]);
	 $px=0;
	 $py=0;
     }
     case "delta1" {
	 $z=0;
	 $pz=1;
	 $x=$rr[$i-1]*cos(2*pi*$rf[$i-1]);
	 $y=$rr[$i-1]*sin(2*pi*$rf[$i-1]);
	 $px=0;
	 $py=0;
     }
     case "delta2" {
	 $z=0;
	 $pz=2;
	 $x=$rr[$i-1]*cos(2*pi*$rf[$i-1]);
	 $y=$rr[$i-1]*sin(2*pi*$rf[$i-1]);
	 $px=0;
	 $py=0;
     }
     case "delta3" {
	 $z=0;
	 $pz=3;
	 $x=$rr[$i-1]*cos(2*pi*$rf[$i-1]);
	 $y=$rr[$i-1]*sin(2*pi*$rf[$i-1]);
	 $px=0;
	 $py=0;
     }
 }
    printf fpw "%f %f %f %f %f %f 1 %d\n",$x,$px,$y,$py,$z,$pz,$i;
}
close(fpw);

$data=[@rr];
my $graph = new GD::Graph::histogram(600,600);
$graph->set( 
x_label         => 'z (sigma)',
y_label         => 'Count',
title           => 'Longitudinal Coordinate Distribution',
x_labels_vertical => 1,
bar_spacing     => 0,
shadow_depth    => 1,
shadowclr       => 'dred',
transparent     => 0,
histogram_bins  => 20,
) or warn $graph->error;
my $gd = $graph->plot($data) or die $graph->error;
        open(IMG, '>histogram.png') or die $!;
        binmode IMG;
        print IMG $gd->png;
