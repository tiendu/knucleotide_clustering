use strict;
use Getopt::Long;

GetOptions(
    "input|i=s" => \my $file_path,
    "knucleotide|k=i" => \my $knu,
    "auto|a=i" => \(my $optimisation = 1),
    "palindromes|p=i" => \(my $use_palindromes = 1),
    "cutoff|c=f" => \my $threshold,
);
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
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
my (%id_sequence_hash, $id);
while (<$input>) {
    chomp;
    if (m/\A>/) {
        ($id) = m/\A>(\S+)/;
    } else {
        $id_sequence_hash{$id} = uc $_;
    };
    last if eof $input;
};
close $input;

our @features;
foreach (kmer_generator($knu)) {
    if ($use_palindromes) {
        push @features, $_ if $_ eq reverse_complement($_);
    } else {
        push @features, $_;
    };
};

my %id_kmer_normalised_count;
for my $id (sort keys %id_sequence_hash) {
    my $no_repeat_sequence = tandem_repeat_remover($id_sequence_hash{$id}, $knu, 3, 0);
    my $kmer_counts = kmer_count_generator($no_repeat_sequence, $knu);
    my $base_frequencies = base_frequency_generator($no_repeat_sequence);
    foreach (@features) {
        my $expected_count = 1;
        my $normalised_usage;
        my $kmer_base_counts = base_count_generator($_);
        for my $base ("A", "T", "G", "C") {
            $expected_count *= ($base_frequencies->{$base} ** $kmer_base_counts->{$base});
        };
        $expected_count *= (length($no_repeat_sequence) - $knu + 1);
        $normalised_usage = $kmer_counts->{$_} / $expected_count;
        $id_kmer_normalised_count{$id}{$_} = $normalised_usage;
    };
};

if ($optimisation == 1) {
    my %silhouette_coefficient = silhouette_coefficient_generator(\%id_kmer_normalised_count);
    my $max = 0;
    my $key;
    print "Threshold\tSilouette coefficient\n";
    foreach (sort keys %silhouette_coefficient) {
        print "$_\t$silhouette_coefficient{$_}\n";
        ($max, $key) = ($silhouette_coefficient{$_}, $_) if $silhouette_coefficient{$_} > $max;
    };
    $threshold = $key;
} elsif ($optimisation == 0) {
    $threshold;
};

my %result = agglomerative_clustering(\%id_kmer_normalised_count, $threshold);
my $i = 1;
foreach (sort keys %result) {
    my @ids = split /,/;
    open my $output, ">:utf8", "cluster_${i}_${file_name}.${file_extension}" or die;
    for my $id (@ids) {
        print "cluster ${i}:\t$id\n";
        print $output ">${id}\n$id_sequence_hash{$id}\n";
    };
    $i++;
    close $output;
};
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
    $sequence = reverse $sequence;
    return $sequence;
};

sub tandem_repeat_remover {
    my ($string, $k, $threshold, $change) = @_;
    foreach (kmer_generator($k)) {
        $string =~ s/($_){$threshold,}/$1 x $change/ge;
    };
    return $string;
};

sub distance_calculator {
    my $distance = 0;
    my %a = %{$_[0]};
    my %b = %{$_[1]};
    my @array = @{$_[2]};
    $distance += ($a{$_} - $b{$_}) ** 2 foreach @array;
    $distance = sqrt($distance);
    return $distance;
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

sub base_count_generator {
    my $sequence = $_[0];
    my @base_list = ("A", "T", "G", "C");
    my $base;
    map {$base->{$_} = 0} @base_list;
    for (my $i = 0; $i <= length($sequence); $i++) {
        $base->{substr($sequence, $i, 1)}++ if exists $base->{substr($sequence, $i, 1)};
    };
    return $base;
};

sub base_frequency_generator {
    my $sequence = $_[0];
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

sub kmer_count_generator {
    my $sequence = $_[0];
    my $k = $_[1];
    my @kmers_list = kmer_generator($k);
    my $kmers;
    map {$kmers->{$_} = 0} @kmers_list;
    for (my $i = 0; $i <= length($sequence) - $k; $i++) {
        $kmers->{substr($sequence, $i, $k)}++ if exists $kmers->{substr($sequence, $i, $k)};
    };
    return $kmers;
};

sub kmer_frequency_generator {
    my $sequence = $_[0];
    my $k = $_[1];
    my @kmers_list = kmer_generator($k);
    my $kmers;
    map {$kmers->{$_} = 0} @kmers_list;
    for (my $i = 0; $i <= length($sequence) - $k; $i++) {
        $kmers->{substr($sequence, $i, $k)}++ if exists $kmers->{substr($sequence, $i, $k)};
    };
    my $total = 0;
    $total += $kmers->{$_} for keys %{$kmers};
    $kmers->{$_} /= $total for keys %{$kmers};
    return $kmers;
};

sub agglomerative_clustering {
    my %data = %{$_[0]};
    my $threshold = $_[1];
    my $size = keys %data;
    my %clusters;
    for (my $i = 1; $i < $size; $i++) {
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
        last if $find->{min} >= $threshold;
        %clusters = %data;
    };
    return %clusters;
};

sub silhouette_coefficient_generator {
    my %data = %{$_[0]};
    my $threshold = 0.1;
    my %parameter_clusters;
    my %parameter_score;
    while ($threshold < sqrt((keys %id_kmer_normalised_count)) / 2) {
        my %result = agglomerative_clustering(\%data, $threshold);
        if (%result) {
            $parameter_clusters{$threshold} = [];
            foreach (sort keys %result) {
                push $parameter_clusters{$threshold}->@*, $_;
            };
            last if keys %result == 2;
        };
        $threshold += 0.1;
    };
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
                            push @intra_distances, distance_calculator($id_kmer_normalised_count{$ids_1[$index_a]}, $id_kmer_normalised_count{$ids_1[$index_b]}, \@::features);
                        };
                        $intra_distance = mean(@intra_distances);
                        $inter_distance = distance_calculator($id_kmer_normalised_count{$ids_1[$index_a]}, $id_kmer_normalised_count{$clusters[$index_2]}, \@::features);
                        push @silhouette_coefficients, ($inter_distance - $intra_distance) / max($inter_distance, $intra_distance);
                    };
                } elsif (scalar @ids_1 == 1 && scalar @ids_2 == 1) {
                    $intra_distance = 0;
                    $inter_distance = distance_calculator($id_kmer_normalised_count{$clusters[$index_1]}, $id_kmer_normalised_count{$clusters[$index_2]}, \@::features);
                    push @silhouette_coefficients, ($inter_distance - $intra_distance) / max($inter_distance, $intra_distance);
                } elsif (scalar @ids_1 == 1 && scalar @ids_2 > 1) {
                    $intra_distance = 0;
                    for my $index_a (0 .. $#ids_2) {
                        push @inter_distances, distance_calculator($id_kmer_normalised_count{$clusters[$index_1]}, $id_kmer_normalised_count{$ids_2[$index_a]}, \@::features);
                    };
                    $inter_distance = mean(@inter_distances);
                    push @silhouette_coefficients, ($inter_distance - $intra_distance) / max($inter_distance, $intra_distance);
                } elsif (scalar @ids_1 > 1 && scalar @ids_2 > 1) {
                    for my $index_a (0 .. $#ids_1) {
                        for my $index_b (1 + $index_a .. $#ids_1) {
                            push @intra_distances, distance_calculator($id_kmer_normalised_count{$ids_1[$index_a]}, $id_kmer_normalised_count{$ids_1[$index_b]}, \@::features);
                        };
                        for my $index_b (0 .. $#ids_2) {
                            push @inter_distances, distance_calculator($id_kmer_normalised_count{$ids_1[$index_a]}, $id_kmer_normalised_count{$ids_2[$index_b]}, \@::features);
                        };
                        $intra_distance = mean(@intra_distances);
                        $inter_distance = mean(@inter_distances);
                        push @silhouette_coefficients, ($inter_distance - $intra_distance) / max($inter_distance, $intra_distance);
                    };
                };
            };
        };
        $parameter_score{$_} = mean(@silhouette_coefficients) unless scalar @silhouette_coefficients == 0;
    };
    return %parameter_score;
};
