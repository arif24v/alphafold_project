#!/usr/bin/perl

print "Ligand_file:\t";
$ligfile = <STDIN>;
chomp $ligfile;

open (FH, $ligfile) || die "Cannot open file $ligfile: $!\n";
@arr_file = <FH>;
close FH;

foreach my $file (@arr_file) {
    chomp $file;
    print "$file\n";
    my @name = split(/\./, $file);
    my $log_file = $file . "_log.log";
    system("\"C:\\Program Files (x86)\\The Scripps Research Institute\\Vina\\vina.exe\" --config conf_vs.txt --ligand \"$file\" --log \"$log_file\"");
}
