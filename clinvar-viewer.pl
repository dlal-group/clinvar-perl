#!/usr/bin/env perl
use warnings;
use List::MoreUtils qw/ uniq /;
use Sort::Fields;

my @file = ();
my $fileToLocate = './variant_summary.clean';
my $fileToLocate2 = './variant_summary.anno';

if (-e $fileToLocate2) {
	print "\n\.anno file is present\n\n";
	open ANNO, "<./variant_summary.anno" or die "$!";
	@file = <ANNO>;
} else {
	if (-e $fileToLocate) {
	    print "\n\.clean file is present\n\n";
	    open CLEAN, "<./variant_summary.clean" or die "$!";
		@file = <CLEAN>;
	} else {
		open CLINVAR, "<./variant_summary.txt" or die "$!";
		open OUT, ">./variant_summary.clean" or die "$!";
	#	open ERROR, ">./variant_summary.error" or die "$!";
		my @clinvar2x = <CLINVAR>;
		@file = grep /\tGRCh37\t/, @clinvar2x;
		unshift(@file, $clinvar2x[0]);
		my @final = ();
		################# clean Clinvar ##################
		print OUT "gene\tvar_type\tconsequence\tsignificance\torigin\ttranscript\taa_ref\tpos_aa\taa_alt\tid\tdiseases\tkey\n";
		foreach my $h (1..$#file) {
			chomp($h);
			my ($id, $type, $mut, $gene, $chr, $ini, $fin, $ref, $alt, $dis, $origin, $sig) = (split /\t/, $file[$h])[0, 1, 2, 4, 18, 19, 20, 21, 22, 13, 14, 6]; # mutacion, gen
			# defitions of extras
			my $consequence = "N/A";
			my $aa_ref = "N/A";
			my $pos_aa = "N/A";
			my $aa_alt = "N/A";
			my $transcript = "N/A";
			my $var = "N/A";
			my $var2 = "N/A";
			# Estructura de la entrada Mut.
			if ($mut =~ /:/) {
				($transcript, $var) = split /:/, $mut, 2;
				$transcript =~ s/\(.*\)//g;
				if ($var =~ /\(*p\./) {
					$var =~ s/.*\(p\.//g;
					$var =~ s/\)//g;
					$var =~ s/\*/Ter/g;
					$aa_ref = $var;
					$aa_ref =~ s/[0-9]+.*//g;
					$var2 = $var;
					$var2 =~ s/[0-9]+/\t/;
					$aa_alt = (split /\t/, $var2)[1];
					$pos_aa = $var;
					$pos_aa =~ s/[a-zA-Z=]+//g;
					if (! defined $aa_ref ) {
						$aa_ref = "N/A";
					} 
					if (! defined $aa_alt) {
						$aa_alt = "N/A";
					} 
					if ($aa_alt eq "=") {
						$aa_alt = $aa_ref;
						$consequence = "synonyomus";
					} elsif ($aa_alt =~ "fs") {
						$consequence = "frameshift";
					} elsif ($aa_alt =~ "del") {
						$consequence = "deletion"; 
					} elsif ($aa_alt ne $aa_ref) {
						$consequence = "missense";
					} 
					if ($aa_alt =~ "Ter") {
						$consequence = "stop_gain";
					}
				} elsif ($var =~ /\Ac\.-/) {
					if ($var =~ /\Ac\.-[0-9]+\+[123][ACTG]/) {
						$consequence = "5-TSS";
					} else {
						$consequence = "5-UTR";				
					}
				} elsif ($var =~ /\Ac\.\*/) {
					if ($var =~ /\Ac\.\*[0-9]+-[123][ACTG]/) {
						$consequence = "3-TTS";
					} else {
						$consequence = "3-UTR";				
					}
				} elsif ($var =~ /\Ac\.[0-9]/) {
					if ($var =~ /\Ac\.[0-9]+\+[12][ACTG]/) {
						$consequence = "splice-donor/aceptor";
					} elsif ($var =~ /\Ac\.[0-9]+-[12][ACTG]/) {
						$consequence = "splice-donor/aceptor";
					} else {
						$consequence = "intronic";				
					}
				} elsif ($var =~ /\Ag\..*del/) {
					$consequence = "deletion";
				} elsif ($var =~ /\Ag\..*dup/) {
					$consequence = "duplication";
				} elsif ($type =~ /\Acopy number/) {
					$consequence = "copy number";
					$pos_aa = $var;
				} elsif ($var =~ /\Ag\./) {
					$consequence = "genomic";
					$pos_aa = $var;
				} elsif ($var =~ /\An\./) {
					$consequence = "non-coding";
				} elsif ($var =~ /\Am\./) {
					$consequence = "mitocondrial";
				} elsif ($type =~ /\Adeletion/) {
					$consequence = "deletion";
					$pos_aa = $var;
				} elsif ($type =~ /\Aduplication/) {
					$consequence = "duplication";
					$pos_aa = $var;
				} elsif ($type =~ /\Acomplex/) {
					$consequence = "complex";
					$pos_aa = $var;
				} else {
	#				print ERROR $type."\t".$var."\n";
				}
		
			} else {
	#			print ERROR "===>".$type."\t".$var."\n";
			}
			my $key = $chr.";".$ini.";".$fin.";".$ref.";".$alt;
		#	print OUT "$key\t$gene\t$type\t$mut\t$consequence\t$sig\t$aa_ref\t$pos_aa\t$aa_alt\t$transcript\t$origin\tCV:$id\t$enf\n"; # all
			print OUT "$gene\t$type\t$consequence\t$sig\t$origin\t$transcript\t$aa_ref\t$pos_aa\t$aa_alt\tCV:$id\t$dis\t$key\n";
		#                0        1        2         3      4          5          6        7        8       9      10    11               
		}
		close CLINVAR;
	}
	close CLEAN;
}
close ANNO;

my @copy = @file;
my $subr = $ARGV[0];
my $argument1 = 0;
my $argument2 = 0;
my @source = ();
my %hash;

if (defined $ARGV[1]) {
	$argument1 = $ARGV[1]; # all keyword
	$argument2 = $ARGV[2]; # query
	if ($argument1 eq "gene") {
		my $header = shift(@copy);
		@file = grep /\A[;,]*$argument2[;,]*\t/, @copy;
		open OUT2, ">./variant_summary.clean.$argument2" or die "$!";
		print OUT2 $header;		
		foreach (@file) {
			print OUT2 $_;
		}
		close OUT2;
	} elsif ($argument1 eq "keyword") {
		my $header = shift(@copy);
		@file = grep /$argument2/, @copy;
		open OUT3, ">./variant_summary.clean.$argument2" or die "$!";
		print OUT3 $header;		
		foreach (@file) {
			print OUT3 $_;
		}
		close OUT3;
	}
}

if (defined $ARGV[0]) {
	if ($subr eq "annotate") {
		open SOURCE, "<./$argument1" or die "$!";
		$argument1 =~ s/\.source//g;
		@source = <SOURCE>;
		# ====================== creando hash
		for my $h (0..$#source) {
			my $fila = shift @source;
			chomp($fila);
			my ($key1, $valor) = (split /\t/, $fila)[0, 1];
			chomp($key1, $valor);
			$hash{$key1} = $valor;
		}
		close SOURCE;
	}
	&$subr;
} else {
	print "\n\tPlease choose one of the following options:
	\n\t1\theader
	2\tsummary
	3\tannotate
	\n";
}

# Header

sub header {
	my @head = split /\t/, shift(@file);
	my $cont = 0;
	foreach my $l (@head) {
		chomp($l);
		$cont++;
		print $cont."\t".$l."\n";
	}
}

# Summary

sub summary {
	# Aqui hay que filtrar la matriz de acuerdo a gen o enfermedad ..........su grep
	# Definicion 
	my @type = ();
	my @consequence = ();
	my @gene = ();
	my @significance = ();
	my @diseases = ();
	# se guardan todas las columnas
	shift(@file);
	foreach my $l (@file) {                                                 # Este clinvar hay que definirlo
		push(@type, (split /\t/, $l)[1]);
		push(@consequence, (split /\t/, $l)[2]);
		my $field0 = (split /\t/, $l)[0];
		my @entrys0 = split /;|, /, $field0;
		foreach my $ci (@entrys0) {
			push(@gene, $ci);
		}
		my $field3 = (split /\t/, $l)[3];
		my @entrys3 = split /;|, /, $field3;
		foreach my $ci (@entrys3) {
			push(@significance, $ci);
		}		
		my $field10 = (split /\t/, $l)[10];
		my @entrys10 = split /;|,/, $field10;
		foreach my $d (@entrys10) {
			push(@diseases, $d);
		}
	}
	# Se extraen los valores unicos de cada columna
	my @u_type = uniq @type;
	my @u_consequence = uniq @consequence;
	my @u_gene = uniq @gene;
	my @u_significance = uniq @significance;
	my @u_diseases = uniq @diseases;

# Aqui se cuenta cuantas entradas tienen la variable
# Type
	my $sum = 0;
	for my $c (0..$#u_type) {
		my $c_type = grep /\A$u_type[$c]\Z/, @type;
		$u_type[$c] = $u_type[$c]."\t".$c_type;
		$sum = $sum + $c_type;
	}
	my @type_sorted = fieldsort '\t', ['-2n'], @u_type;
	push(@type_sorted, "TOTAL\t$sum"); 
# Consequence
	my $sum2 = 0;
	for my $c (0..$#u_consequence) {
		my $c_consequence = grep /\A$u_consequence[$c]\Z/, @consequence;
		$u_consequence[$c] = $u_consequence[$c]."\t".$c_consequence;
		$sum2 = $sum2 + $c_consequence;
	}
	my @consequence_sorted = fieldsort '\t', ['-2n'], @u_consequence;
	push(@consequence_sorted, "TOTAL\t$sum2"); 
# Significance
	my $sum3 = 0;
	for my $c (0..$#u_significance) {
		my $c_significance = grep /\A$u_significance[$c]\Z/, @significance;
		$u_significance[$c] = $u_significance[$c]."\t".$c_significance;
		$sum3 = $sum3 + $c_significance;
	}
	my @sig_sorted = fieldsort '\t', ['-2n'], @u_significance;
	push(@sig_sorted, "TOTAL\t$sum3"); 
# Se imprime
	print "\n <=== Summary for $argument1: $argument2 ===>\n"; # aqui colocar el nombre del query
	print "\n### TYPES = ".($#u_type+1)."\n";
	foreach (@type_sorted) {
		print "\t".$_."\n";
	}
	print "\n### CONSEQUENCE = ".($#u_consequence+1)."\n";
	foreach (@consequence_sorted) {
		print "\t".$_."\n";
	}
	print "\n### CLINICAL SIGNIFICANCE = ".($#u_significance+1)."\n";
	foreach (@sig_sorted) {
		print "\t".$_."\n";
	}
######### Only if filter provided
	if ($argument1 ne "0") {
#	# Genes
		my $sum4 = 0;
		for my $c (0..$#u_gene) {
			my $c_gene = grep /\A$u_gene[$c]\Z/, @gene;
			$u_gene[$c] = $u_gene[$c]."\t".$c_gene;
			$sum4 = $sum4 + $c_gene;
		}
		my @gene_sorted = fieldsort '\t', ['-2n'], @u_gene;
		push(@gene_sorted, "TOTAL\t$sum4"); 
		print "\n### GENES = ".($#u_gene+1)."\n";
		for my $g (0..$#gene_sorted) {
			print "\t".$gene_sorted[$g]."\n";
		}
		print "\t...\n";
		print "\tTOTAL\t$sum4\n";
#	 Diseases
		my $sum5 = 0;
		for my $c (0..$#u_diseases) {
			my $c_diseases = grep /\A$u_diseases[$c]\Z/, @diseases;
			$u_diseases[$c] = $u_diseases[$c]."\t".$c_diseases;
			$sum5 = $sum5 + $c_diseases;
		}
		my @diseases_sorted = fieldsort '\t', ['-2n'], @u_diseases;
		push(@diseases_sorted, "TOTAL\t$sum5"); 
		print "\n### DISEASES = ".($#u_diseases+1)."\n";
		for my $g (0..$#diseases_sorted) {
			print "\t".$diseases_sorted[$g]."\n";
		}
		print "\t...\n";
		print "\tTOTAL\t$sum5\n";
	}
}

# Annotate

sub annotate {
	my $header = shift(@file);
	chomp($header);
	my $annoheader = $header."\t$argument1\n";
	open OUT4, ">./variant_summary.anno" or die "$!";
	print OUT4 $annoheader;
	for my $l (0..$#file) {
		chomp($file[$l]);		
		my $key = (split /\t/, $file[$l])[11];
		if (exists $hash{$key}) {
			print OUT4 $file[$l]."\t".$hash{$key}."\n";
		} else {
			print OUT4 $file[$l]."\tN/A\n";
		}
	}
	print "\nannotation successfull\n";
}