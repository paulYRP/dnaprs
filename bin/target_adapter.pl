#!/usr/bin/env perl
use strict;
use warnings;
use Text::ParseWords qw(parse_line);

sub usage {
    die "Usage: target_adapter.pl finalreport <report> <manifest.csv> <output-prefix>\n" .
        "   or: target_adapter.pl annotate-pvar <input.pvar> <manifest.csv|''> <marker-map.tsv|''> <output.pvar> <decisions.tsv>\n";
}

sub complement {
    my ($value) = @_;
    $value //= '';
    $value =~ tr/ACGT/TGCA/;
    return $value;
}

sub read_manifest {
    my ($path) = @_;
    my %result;
    return \%result if !defined($path) || $path eq '';
    open my $handle, '<', $path or die "Cannot read assay manifest $path: $!\n";
    my @header;
    my %column;
    my $in_assay = 0;
    while (my $line = <$handle>) {
        $line =~ s/\r?\n$//;
        if (!$in_assay) {
            $in_assay = 1 if $line eq '[Assay]';
            next;
        }
        if (!@header) {
            @header = parse_line(',', 0, $line);
            @column{@header} = (0 .. $#header);
            for my $required (qw(Name IlmnStrand SNP Chr MapInfo RefStrand)) {
                die "Assay manifest is missing $required\n" if !exists $column{$required};
            }
            next;
        }
        next if $line eq '' || $line =~ /^\[/;
        my @value = parse_line(',', 0, $line);
        my $name = $value[$column{Name}] // '';
        next if $name eq '' || exists $result{$name};
        my $chromosome = $value[$column{Chr}] // '';
        my $position = $value[$column{MapInfo}] // '';
        my $strand = uc($value[$column{IlmnStrand}] // '');
        my $ref_strand = $value[$column{RefStrand}] // '';
        my $snp = uc($value[$column{SNP}] // '');
        my ($allele_a, $allele_b) = $snp =~ /^\[([ACGT])\/([ACGT])\]$/;
        # Convert the assay-coded alleles to TOP first, then TOP to the
        # reference-forward strand used by PLINK and the GRCh build.
        if (defined $allele_a && $strand eq 'BOT') {
            $allele_a = complement($allele_a);
            $allele_b = complement($allele_b);
        }
        if (defined $allele_a && $ref_strand eq '-') {
            $allele_a = complement($allele_a);
            $allele_b = complement($allele_b);
        }
        $result{$name} = {
            chromosome => $chromosome,
            position => $position,
            allele_a => $allele_a // '',
            allele_b => $allele_b // '',
            complement_calls => $ref_strand eq '-' ? 1 : 0,
        };
    }
    close $handle;
    return \%result;
}

sub read_marker_map {
    my ($path) = @_;
    my %result;
    return \%result if !defined($path) || $path eq '';
    open my $handle, '<', $path or die "Cannot read marker map $path: $!\n";
    my $header = <$handle> // die "Marker map is empty: $path\n";
    $header =~ s/\r?\n$//;
    my @header = split /\t/, $header, -1;
    my %column;
    @column{@header} = (0 .. $#header);
    die "Marker map requires source_id and new_id columns\n"
        if !exists $column{source_id} || !exists $column{new_id};
    while (my $line = <$handle>) {
        $line =~ s/\r?\n$//;
        next if $line eq '';
        my @value = split /\t/, $line, -1;
        my $source = $value[$column{source_id}] // '';
        next if $source eq '' || exists $result{$source};
        my %record;
        for my $field (qw(new_id chr pos ref alt)) {
            $record{$field} = exists $column{$field} ? ($value[$column{$field}] // '') : '';
        }
        $result{$source} = \%record;
    }
    close $handle;
    return \%result;
}

sub finalreport_to_ped {
    my ($report, $manifest_path, $prefix) = @_;
    my $manifest = read_manifest($manifest_path);
    die "No assay records were read from $manifest_path\n" if !%{$manifest};

    open my $input, '<', $report or die "Cannot read GenomeStudio report $report: $!\n";
    open my $ped, '>', "$prefix.ped" or die "Cannot write $prefix.ped: $!\n";
    open my $map, '>', "$prefix.map" or die "Cannot write $prefix.map: $!\n";
    my %column;
    my $in_data = 0;
    my $sample = '';
    my $marker_count = 0;
    my $expected_count;
    while (my $line = <$input>) {
        $line =~ s/\r?\n$//;
        if (!$in_data) {
            $in_data = 1 if $line eq '[Data]';
            next;
        }
        if (!%column) {
            my @header = split /\t/, $line, -1;
            @column{@header} = (0 .. $#header);
            for my $required ('SNP Name', 'Sample ID', 'Allele1 - Top', 'Allele2 - Top') {
                die "GenomeStudio report is missing $required\n" if !exists $column{$required};
            }
            next;
        }
        next if $line eq '';
        my @value = split /\t/, $line, -1;
        my $marker = $value[$column{'SNP Name'}] // '';
        my $current_sample = $value[$column{'Sample ID'}] // '';
        next if !exists $manifest->{$marker};
        my $record = $manifest->{$marker};
        next if $record->{chromosome} !~ /^(?:[1-9]|1[0-9]|2[0-2]|X|Y|XY|MT)$/i;
        next if $record->{position} !~ /^\d+$/ || $record->{position} < 1;
        next if $record->{allele_a} eq '' || $record->{allele_b} eq '';

        if ($sample ne $current_sample) {
            if ($sample ne '') {
                print {$ped} "\n";
                $expected_count //= $marker_count;
                die "GenomeStudio samples do not share one marker order/count\n" if $marker_count != $expected_count;
            }
            $sample = $current_sample;
            $marker_count = 0;
            print {$ped} join(' ', $sample, $sample, 0, 0, 0, -9);
        }
        my $allele1 = uc($value[$column{'Allele1 - Top'}] // '0');
        my $allele2 = uc($value[$column{'Allele2 - Top'}] // '0');
        if ($record->{complement_calls}) {
            $allele1 = complement($allele1);
            $allele2 = complement($allele2);
        }
        $allele1 = '0' if $allele1 !~ /^[ACGT]$/;
        $allele2 = '0' if $allele2 !~ /^[ACGT]$/;
        print {$ped} " $allele1 $allele2";
        if (!defined $expected_count) {
            my %chromosome_code = (X => 23, Y => 24, XY => 25, MT => 26);
            my $chromosome = uc($record->{chromosome});
            $chromosome = $chromosome_code{$chromosome} if exists $chromosome_code{$chromosome};
            print {$map} join("\t", $chromosome, $marker, 0, $record->{position}), "\n";
        }
        $marker_count++;
    }
    if ($sample ne '') {
        print {$ped} "\n";
        $expected_count //= $marker_count;
        die "GenomeStudio samples do not share one marker count\n" if $marker_count != $expected_count;
    }
    close $input;
    close $ped;
    close $map;
    die "GenomeStudio conversion produced no samples or markers\n" if !defined($expected_count) || $expected_count == 0;
}

sub annotate_pvar {
    my ($input_path, $manifest_path, $map_path, $output_path, $decision_path) = @_;
    my $manifest = read_manifest($manifest_path);
    my $marker_map = read_marker_map($map_path);
    open my $input, '<', $input_path or die "Cannot read $input_path: $!\n";
    open my $output, '>', $output_path or die "Cannot write $output_path: $!\n";
    open my $decision, '>', $decision_path or die "Cannot write $decision_path: $!\n";
    print {$decision} join("\t", qw(source_id final_id source_chr source_pos final_chr final_pos source_ref source_alt final_ref final_alt decision reason)), "\n";
    my %seen_id;
    while (my $line = <$input>) {
        if ($line =~ /^##/) {
            print {$output} $line;
            next;
        }
        $line =~ s/\r?\n$//;
        if ($line =~ /^#CHROM\t/) {
            print {$output} "$line\n";
            next;
        }
        my @value = split /\t/, $line, -1;
        die "Unexpected PVAR row: $line\n" if @value < 5;
        my ($source_chr, $source_pos, $source_id, $source_ref, $source_alt) = @value[0 .. 4];
        my ($final_chr, $final_pos, $final_id, $final_ref, $final_alt) =
            ($source_chr, $source_pos, $source_id, uc($source_ref), uc($source_alt));
        my (@reason, $matched);
        if (exists $marker_map->{$source_id}) {
            my $record = $marker_map->{$source_id};
            $final_id = $record->{new_id} if $record->{new_id} ne '';
            $final_chr = $record->{chr} if $record->{chr} ne '';
            $final_pos = $record->{pos} if $record->{pos} ne '';
            $final_ref = uc($record->{ref}) if $record->{ref} ne '';
            $final_alt = uc($record->{alt}) if $record->{alt} ne '';
            push @reason, 'explicit marker map';
            $matched = 1;
        } elsif (exists $manifest->{$source_id}) {
            my $record = $manifest->{$source_id};
            $final_chr = $record->{chromosome} if $record->{chromosome} ne '';
            $final_pos = $record->{position} if $record->{position} ne '';
            my @assay = ($record->{allele_a}, $record->{allele_b});
            if (($final_ref eq '' || $final_ref eq '.' || $final_ref eq '0') && $final_alt =~ /^[ACGT]$/) {
                my ($other) = grep { $_ ne $final_alt } @assay;
                $final_ref = $other if defined $other;
            }
            if (($final_alt eq '' || $final_alt eq '.' || $final_alt eq '0') && $final_ref =~ /^[ACGT]$/) {
                my ($other) = grep { $_ ne $final_ref } @assay;
                $final_alt = $other if defined $other;
            }
            push @reason, 'Illumina assay manifest';
            $matched = 1;
        }
        my $decision_value = $matched ? 'CORRECTED_OR_VERIFIED' : 'REVIEW';
        push @reason, 'no marker annotation supplied' if !$matched;
        if ($seen_id{$final_id}++) {
            $decision_value = 'REVIEW';
            push @reason, 'duplicate final identifier';
        }
        @value[0 .. 4] = ($final_chr, $final_pos, $final_id, $final_ref, $final_alt);
        print {$output} join("\t", @value), "\n";
        print {$decision} join("\t", $source_id, $final_id, $source_chr, $source_pos, $final_chr,
            $final_pos, $source_ref, $source_alt, $final_ref, $final_alt, $decision_value,
            join('; ', @reason)), "\n";
    }
    close $input;
    close $output;
    close $decision;
}

my $action = shift @ARGV // usage();
if ($action eq 'finalreport') {
    @ARGV == 3 or usage();
    finalreport_to_ped(@ARGV);
} elsif ($action eq 'annotate-pvar') {
    @ARGV == 5 or usage();
    annotate_pvar(@ARGV);
} else {
    usage();
}
