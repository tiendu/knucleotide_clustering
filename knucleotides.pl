use strict;
use warnings;
use Getopt::Long;

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package IdSequence;

sub new {
    my ($class, $id, $sequence) = @_;
    my $self = {
        _id => $id,
        _sequence => $sequence,
    };
    return bless($self, $class);
};

sub base_frequency {
    my ($self) = @_;
    my $sequence = $self->{_sequence};
    my @base_list = ("A", "T", "G", "C");
    my $base;
    map {$base->{$_} = 0} @base_list;
    for (my $i = 0; $i <= length($sequence); $i++) {
        $base->{substr($sequence, $i, 1)}++ if exists $base->{substr($sequence, $i, 1)};
    };
    my $total = 0;
    $total += $base->{$_} for keys %{$base};
    $base->{$_} /= $total for keys %{$base};
    return $base;
};

sub kmer_frequency {
    my ($self, $k) = @_;
    my $sequence = $self->{_sequence};
    my @bases_1 = ("A", "T", "G", "C");
    my @bases_2 = @bases_1;
    for (my $i = 1; $i < $k; $i++) {
        my @temporary;
        for my $base_1 (@bases_1) {
            for my $base_2 (@bases_2) {
                push @temporary, "$base_1" . "$base_2";
            };
        };
        undef @bases_2;
        @bases_2 = @temporary;
    };
    my $kmers;
    map {$kmers->{$_} = 0} @bases_2;
    for (my $i = 0; $i <= length($sequence) - $k; $i++) {
        $kmers->{substr($sequence, $i, $k)}++ if exists $kmers->{substr($sequence, $i, $k)};
    };
    my $total = 0;
    $total += $kmers->{$_} for keys %{$kmers};
    $kmers->{$_} /= $total for keys %{$kmers};
    return $kmers;
};

sub normalised_kmer_frequency {
    my ($self, $k) = @_;
    my $base_freq = $self->base_frequency();
    my $kmer_freq = $self->kmer_frequency($k);
    my %normal_freq;
    for my $knu (keys %{$kmer_freq}) {
        $normal_freq{$knu} = $kmer_freq->{$knu};
        my @nus = split //, $knu;
        foreach (@nus) {
            my $test = $base_freq->{$_};
            if ($base_freq->{$_} != 0) {
                $normal_freq{$knu} /= $base_freq->{$_};
            };
        };
    };
    return %normal_freq;
};
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package main;
GetOptions(
    "input|i=s" => \my $file_path,
    "knucleotide|k=i" => \my $knu,
    "auto|a=i" => \(my $optimisation = 1),
    "palindromes|p=i" => \(my $use_palindromes = 1),
    "cutoff|c=f" => \my $threshold,
);

if (!$file_path) {
    print "$0: Input required!\n";
    exit 1;
} elsif (!$knu) {
    print "$0: k-nucleotide required!\n";
    exit 1;
} elsif ($optimisation == 0 && !$threshold) {
    print "$0: Manual clustering, threshold required!\n";
    exit 1;
};

my $file_path_cp = $file_path;
$file_path_cp =~ s/.*\///;
my ($file_name, $file_extension) = $file_path_cp =~ /^(.+)\.([^.]+)$/;

print "Input: ${file_path}\n";
if ($optimisation == 1) {
    print "Parameters:\n\t- k-Nucleotide: ${knu}\n\t- Silhouette optimisation: ${optimisation}\n\t- Use palindromes: ${use_palindromes}\n\t- Dendrogram cutoff: Auto\n";
} else {
    print "Parameters:\n\t- k-Nucleotide: ${knu}\n\t- Silhouette optimisation: ${optimisation}\n\t- Use palindromes: ${use_palindromes}\n\t- Dendrogram cutoff: ${threshold}\n";
}
print "=" x 30 . "\n";
print "Running now! Please wait...\n";

open my $input, "<:utf8", $file_path or die;
my (@IdSequence_array, $id);
while (<$input>) {
    chomp;
    if (m/\A>/) {
        ($id) = m/\A>(\S+)/;
    } else {
        my $a_IdSequence = IdSequence->new($id, uc $_);
        push @IdSequence_array, $a_IdSequence;
    };
    last if eof $input;
};
close $input;

my %id_kmer_normalised_frequency;
for my $IdSequence (@IdSequence_array) {
    $IdSequence->{_sequence} = tandem_repeat_remover($IdSequence->{_sequence}, $knu, 3, 0);
    my %normalised_frequency = $IdSequence->normalised_kmer_frequency($knu);
    foreach (keys %normalised_frequency) {
        $id_kmer_normalised_frequency{$IdSequence->{_id}}{$_} = $normalised_frequency{$_};
    };
};

our @features;
foreach (kmer_generator($knu)) {
    if ($use_palindromes) {
        push @features, $_ if $_ eq reverse_complement($_);
    } else {
        push @features, $_;
    };
};

if ($optimisation == 1) {
    my %silhouette_coefficient = silhouette_coefficient_generator(\%id_kmer_normalised_frequency);
    my $max = 0;
    my $key;
    open my $report, ">:utf8", "silhouette_coef_${file_name}.txt" or die;
    print $report "Threshold\tNo. of clusters\tSilouette coefficient\n";
    for my $parameter (sort keys %silhouette_coefficient) {
        my $number_of_clusters = $silhouette_coefficient{$parameter}[0];
        my $coefficient = $silhouette_coefficient{$parameter}[1];
        print $report "$parameter\t$number_of_clusters\t$coefficient\n";
        ($max, $key) = ($coefficient, $parameter) if $coefficient > $max;
    };
    $threshold = $key;
    close $report;
} elsif ($optimisation == 0) {
    $threshold = $threshold;
};

my %result = agglomerative_clustering(\%id_kmer_normalised_frequency, $threshold);
my $i = 1;
open my $output_1, ">:utf8", "clusters_${file_name}.txt" or die;
foreach (sort keys %result) {
    my @ids = split /,/;
    open my $output_2, ">:utf8", "cluster_${i}_${file_name}.${file_extension}" or die;
    print $output_1 "#" x 30 . "\n";
    for my $id (@ids) {
        print $output_1 "cluster ${i}:\t$id\n";
        for my $IdSequence (@IdSequence_array) {
            print $output_2 ">${id}" . "\n" . $IdSequence->{_sequence} . "\n" if $id eq $IdSequence->{_id};
        };
    };
    $i += 1;
    close $output_2;
};
close $output_1;
print "=" x 30 . "\n";
print "Clustering done!\n";
print "Output:\n" . "- Optimisation report: silhouette_coef_${file_name}.txt\n" . "- Cluster information: clusters_${file_name}.txt\n" . "- Files contain the sequences in FASTA format: cluster_n_${file_name}.${file_extension}\n" if $optimisation == 1;
print "Output:\n" . "- Cluster information: clusters_${file_name}.txt\n" . "- Files contain the sequences in FASTA format: cluster_n_${file_name}.${file_extension}\n" if $optimisation == 0;
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
sub mean {
    my @data = @_;
    my $sum;
    $sum += $_ for @data;    
    return $sum / @data;
};

sub max {
    my $max = shift;
    foreach (@_) {
        $max = $_ if $_ > $max;
    };
    return $max;
}

sub min {
    my $min = shift;
    foreach (@_) {
        $min = $_ if $_ < $min;
    };
    return $min;
};

sub palindrome {
    my $word = $_;
    my $reversed_word = reverse $word;
    if ($word eq $reversed_word) {
        return 1;
    } else {
        return 0;
    };
};

sub reverse_complement {
    my $sequence = $_[0];
    $sequence =~ tr/ATGC/TACG/;
    return reverse($sequence);
};

sub tandem_repeat_remover {
    my ($string, $k, $threshold, $change) = @_;
    foreach (kmer_generator($k)) {
        $string =~ s/($_){$threshold,}/$1 x $change/ge;
    };
    return $string;
};

sub kmer_generator {
    my ($k) = @_;
    my @bases_1 = ("A", "T", "G", "C");
    my @bases_2 = @bases_1;
    for (my $i = 1; $i < $k; $i++) {
        my @temporary;
        for my $base_1 (@bases_1) {
            for my $base_2 (@bases_2) {
                push @temporary, "$base_1" . "$base_2";
            };
        };
        undef @bases_2;
        @bases_2 = @temporary;
    };
    return @bases_2;
};

sub distance_calculator {
    my $distance = 0;
    my %a = %{$_[0]};
    my %b = %{$_[1]};
    my @array = @{$_[2]};
    $distance += ($a{$_} - $b{$_}) ** 2 foreach @array;
    return sqrt($distance);
};

sub agglomerative_clustering {
    my %data = %{$_[0]};
    my $threshold = $_[1];
    my %clusters;
    for (my $i = 1; $i < (keys %data); $i++) {
        my (%distances, $find, @keys);
        @keys = sort keys %data;
        for my $index_1 (0 .. $#keys) {
            for my $index_2 (1 + $index_1 .. $#keys) {
                my ($distance, $key_1, $key_2) = (0, $keys[$index_1], $keys[$index_2]);
                $distance = distance_calculator($data{$key_1}, $data{$key_2}, \@::features);
                $distances{$key_1}{$key_2} = $distance;
                $find->{min} = $distance unless $find->{min};
                $find->{key} = [$key_1, $key_2] unless $find->{key};
                if ($find->{min} > $distance) {
                    $find->{min} = $distance;
                    $find->{key} = [$key_1, $key_2];
                };
            };
        };
        my ($key_1, $key_2) = $find->{key}->@*;
        $data{"$key_1,$key_2"}{$_} = mean($data{$key_1}{$_}, $data{$key_2}{$_}) foreach @::features;
        delete @data{($key_1, $key_2)};
        %clusters = %data;
        last if $find->{min} >= $threshold;
    };
    return %clusters;
};

sub silhouette_coefficient_generator {
    my %data = %{$_[0]};
    my $threshold = 0.1;
    my %parameter_clusters;
    while (1) {
        my %result = agglomerative_clustering(\%data, $threshold);
        if (%result) {
            $parameter_clusters{$threshold} = [];
            foreach (sort keys %result) {
                push $parameter_clusters{$threshold}->@*, $_;
            };
            last if keys %result == (keys %data) / 2;
        };
        $threshold += 0.1;
    };
    my %parameter_score;
    foreach (sort keys %parameter_clusters) {
        my @silhouette_coefficients;
        my @clusters = $parameter_clusters{$_}->@*;
        for my $index_1 (0 .. $#clusters) {
            for my $index_2 (1 + $index_1 .. $#clusters) {
                my @ids_1 = split /,/, $clusters[$index_1];
                my @ids_2 = split /,/, $clusters[$index_2];
                my ($intra_distance, $inter_distance, @intra_distances, @inter_distances);
                if (scalar @ids_1 > 1 && scalar @ids_2 == 1) {
                    for my $index_a (0 .. $#ids_1) {
                        for my $index_b (1 + $index_a .. $#ids_1) {
                            push @intra_distances, distance_calculator($data{$ids_1[$index_a]}, $data{$ids_1[$index_b]}, \@::features);
                        };
                        $intra_distance = mean(@intra_distances);
                        $inter_distance = distance_calculator($data{$ids_1[$index_a]}, $data{$clusters[$index_2]}, \@::features);
                    };
                } elsif (scalar @ids_1 == 1 && scalar @ids_2 == 1) {
                    $intra_distance = 0;
                    $inter_distance = distance_calculator($data{$clusters[$index_1]}, $data{$clusters[$index_2]}, \@::features);
                } elsif (scalar @ids_1 == 1 && scalar @ids_2 > 1) {
                    $intra_distance = 0;
                    for my $index_a (0 .. $#ids_2) {
                        push @inter_distances, distance_calculator($data{$clusters[$index_1]}, $data{$ids_2[$index_a]}, \@::features);
                    };
                    $inter_distance = mean(@inter_distances);
                } elsif (scalar @ids_1 > 1 && scalar @ids_2 > 1) {
                    for my $index_a (0 .. $#ids_1) {
                        for my $index_b (1 + $index_a .. $#ids_1) {
                            push @intra_distances, distance_calculator($data{$ids_1[$index_a]}, $data{$ids_1[$index_b]}, \@::features);
                        };
                        for my $index_b (0 .. $#ids_2) {
                            push @inter_distances, distance_calculator($data{$ids_1[$index_a]}, $data{$ids_2[$index_b]}, \@::features);
                        };
                        $intra_distance = mean(@intra_distances);
                        $inter_distance = mean(@inter_distances);
                    };
                };
                push @silhouette_coefficients, ($inter_distance - $intra_distance) / max($inter_distance, $intra_distance);
            };
        };
        $parameter_score{$_} = [scalar @clusters, mean(@silhouette_coefficients)] unless scalar @silhouette_coefficients == 0;
    };
    return %parameter_score;
};
