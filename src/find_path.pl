#!usr/bin/perl      -w

use strict;
use Data::Dumper;
use Storable qw(lock_nstore);
use Storable qw(lock_retrieve);

############################################################################################################################
# global variable
my $max_level=31;
my $dir="./GEM/";
if (@ARGV > 0) {
    $dir=$ARGV[0];
}

my $num=0;
my $stop=0;
my %medium_end=();#%{lock_retrieve("gem_info/storable/end_reactant/saccharides_and_amino_acids.storable")};
$medium_end{'H'}=1;
$medium_end{'H2O'}=1;
$medium_end{'H2O H2O'}=1;
$medium_end{'UTP'}=1;
$medium_end{'UTP C9H11N2O15P3'}=1;
$medium_end{'ATP C10H12N5O13P3'}=1;
$medium_end{'ATP'}=1;
$medium_end{'Phosphate'}=1;
$medium_end{'Nicotinamide adenine dinucleotide phosphate - reduced'}=1;
my %add_reaction=(); # added reactions
#$add_reaction{'Sucrose 6-phosphate'}="D-Glucose";
#$add_reaction{'Malonate semialdehyde'}="malonyl-CoA";  
#$add_reaction{'malonyl-CoA'}="D-Fructose 6-phosphate";  
#$add_reaction{'D-Fructose 6-phosphate'}="D-Fructose";  
#$add_reaction{'D-Fructose'}="Sucrose C12H22O11";  
#$add_reaction{'Sucrose C12H22O11'}="D-Fructose";  
#$add_reaction{'Sucrose C12H22O11'}="D-Glucose";  
#$add_reaction{'Malonate semialdehyde'}="Sucrose C12H22O11";  
#my %metabolite=%{lock_retrieve("specific_bacteria_gem_info/storable/metabolite.storable")};
#$metabolite{'M_sucr_e'}{'name'}="Sucrose C12H22O11";
#$metabolite{'M_sucr_p'}{'name'}="Sucrose C12H22O11";
#$metabolite{'M_sucr_c'}{'name'}="Sucrose C12H22O11";
my %basic_molecules; # currency metabolites and other basic molecules
$basic_molecules{'ADP C10H12N5O10P2'}=1;
$basic_molecules{'ADP'}=1;
$basic_molecules{'H+'}=1;
$basic_molecules{'Phosphate'}=1;
$basic_molecules{'ATP C10H12N5O13P3'}=1;
$basic_molecules{'ATP'}=1;
$basic_molecules{'UTP C9H11N2O15P3'}=1;
$basic_molecules{'UTP'}=1;
$basic_molecules{'H2O H2O'}=1;
$basic_molecules{'H2O'}=1;
$basic_molecules{'Nicotinamide adenine dinucleotide - reduced'}=1;
$basic_molecules{'Nicotinamide adenine dinucleotide phosphate'}=1;
$basic_molecules{'Nicotinamide adenine dinucleotide phosphate - reduced'}=1;
$basic_molecules{'Nicotinamide adenine dinucleotide'}=1;
$basic_molecules{'O2 O2'}=1;
$basic_molecules{'O2'}=1;
$basic_molecules{'Hydrogen sulfide'}=1;
$basic_molecules{'Coenzyme A'}=1;
# $basic_molecules{'Acetyl-CoA'}=1;
$basic_molecules{'Diphosphate'}=1;
$basic_molecules{'AMP C10H12N5O7P'}=1;
$basic_molecules{'AMP'}=1;
$basic_molecules{'Ammonium'}=1;
$basic_molecules{'Co2+'}=1;
$basic_molecules{'Potassium'}=1;
$basic_molecules{'Potassium'}=1;
$basic_molecules{'1-(5-Phosphoribosyl)-ATP'}=1;
$basic_molecules{'1-alkyl 2-lysoglycerol 3-phosphocholine'}=1;
$basic_molecules{'UDP C9H11N2O12P2'}=1;
$basic_molecules{'UDP'}=1;
$basic_molecules{'CMP C9H12N3O8P'}=1;
$basic_molecules{'CMP'}=1;
$basic_molecules{'CTP C9H12N3O14P3'}=1;
$basic_molecules{'CTP'}=1;
$basic_molecules{'GTP C10H12N5O14P3'}=1;
$basic_molecules{'GTP'}=1;
$basic_molecules{'DADP C10H12N5O9P2'}=1;
$basic_molecules{'DADP'}=1;
$basic_molecules{'DATP C10H12N5O12P3'}=1;
$basic_molecules{'DATP'}=1;
$basic_molecules{'DTDP C10H13N2O11P2'}=1;
$basic_molecules{'DTDP'}=1;
$basic_molecules{'Flavin adenine dinucleotide oxidized'}=1;
$basic_molecules{'Flavin adenine dinucleotide'}=1;
$basic_molecules{'IDP C10H11N4O11P2'}=1;
$basic_molecules{'IDP'}=1;
$basic_molecules{'Sodium'}=1;
$basic_molecules{'DADP C10H12N5O9P2'}=1;
$basic_molecules{'DADP'}=1;
$basic_molecules{'ITP C10H11N4O14P3'}=1;
$basic_molecules{'ITP'}=1;
$basic_molecules{'IMP C10H11N4O8P'}=1;
$basic_molecules{'IMP'}=1;
$basic_molecules{'Biotin'}=1;
$basic_molecules{'DGMP C10H12N5O7P'}=1;
$basic_molecules{'DGMP'}=1;
$basic_molecules{'DTDP C10H13N2O11P2'}=1;
$basic_molecules{'DTDP'}=1;
$basic_molecules{'DAMP C10H12N5O6P'}=1;
$basic_molecules{'DAMP'}=1;
$basic_molecules{'DCMP C9H12N3O7P'}=1;
$basic_molecules{'DCMP'}=1;
$basic_molecules{'DATP C10H12N5O12P3'}=1;
$basic_molecules{'DATP'}=1;
$basic_molecules{'DTMP C10H13N2O8P'}=1;
$basic_molecules{'DTMP'}=1;
$basic_molecules{'UMP C9H11N2O9P'}=1;
$basic_molecules{'UMP'}=1;
$basic_molecules{'Fe2+ mitochondria'}=1;
$basic_molecules{'DGTP C10H12N5O13P3'}=1;
$basic_molecules{'DGTP'}=1;
$basic_molecules{'DGDP C10H12N5O10P2'}=1;
$basic_molecules{'DGDP'}=1;
$basic_molecules{'DTTP C10H13N2O14P3'}=1;
$basic_molecules{'DTTP'}=1;
$basic_molecules{'Igg lc[r]'}=1;
$basic_molecules{'Igg hc[r]'}=1;
$basic_molecules{'GDP C10H12N5O11P2'}=1;
$basic_molecules{'GDP'}=1;
$basic_molecules{'EpoIm[r]'}=1;
$basic_molecules{'Igg[g]'}=1;
$basic_molecules{'I c'}=1;
$basic_molecules{'Epo[g]'}=1;
$basic_molecules{'Peptide[c]'}=1;
$basic_molecules{'Sodium'}=1;
$basic_molecules{'Nickle'}=1;
$basic_molecules{'H+'}=1;
$basic_molecules{'Copper'}=1;
$basic_molecules{'Calcium'}=1;
$basic_molecules{'Zinc'}=1;
$basic_molecules{'Sulfate'}=1;
$basic_molecules{'Fe3+'}=1;
$basic_molecules{'Magnesium'}=1;
$basic_molecules{'Manganese'}=1;
$basic_molecules{'Nickel'}=1;
$basic_molecules{'Photon'}=1;
$basic_molecules{'Silver'}=1;
$basic_molecules{'Rifampin'}=1;
$basic_molecules{'Iron'}=1;
$basic_molecules{'GMP C10H12N5O8P'}=1;
$basic_molecules{'GMP'}=1;
$basic_molecules{'Cu+'}=1;
$basic_molecules{'Cadmium'}=1;
$basic_molecules{'Zinc'}=1;
$basic_molecules{'DCTP C9H12N3O13P3'}=1;
$basic_molecules{'DCTP'}=1;
$basic_molecules{'DUTP C9H11N2O14P3'}=1;
$basic_molecules{'DUTP'}=1;
$basic_molecules{'DUMP C9H11N2O8P'}=1;
$basic_molecules{'DUMP'}=1;
#$basic_molecules{'CO2 CO2'}=1;
#$basic_molecules{'CO2'}=1;
$basic_molecules{''}=1;
my %target_media;
#$target_media{'CO2 CO2'}=1;
#$target_media{'D-Glucose'}=1;
#$target_media{'Sucrose C12H22O11'}=1;

my %bigg_not_bacteria; #exclude non-bacteria GEMs in BiGG database (deprecated)
$bigg_not_bacteria{'RECON1'}=1;
$bigg_not_bacteria{'Recon3D'}=1;
$bigg_not_bacteria{'iAT_PLT_636'}=1;
$bigg_not_bacteria{'iCHOv1'}=1;
$bigg_not_bacteria{'iMM1415'}=1;
my %new_carve_me_not_bacteria;
my @csources;
############################################################################################################################
## main program

new_carve_me_read_and_analyze_specific_bacteria_gem(); # analyze xml GEM file
new_carve_me_find_path_from_secretion_to_intake_for_each_bacteria(); # find pathways
#new_carve_me_analyze_path_result("3-Hydroxypropanoate","Sucrose","CO2"); # find combinations for specific pathways (deprecated)

#new_carve_me_deep_analysis_of_bacteria_internal_pathways_and_link_target_secretion_to_target_intake("Formate","L-Methionine","ASF356"); # find specific pathways of certein intake and secretion

#my %bac_reactions=%{lock_retrieve("/data/ylliao/ebiota/bigg/reaction_info_for_every_bacteria/storable/iJB785.storable")}; # check reaction of specific GEM
#print_reaction("R_RBPCcx",\%bac_reactions);

############################################################################################################################
## functions
sub new_carve_me_read_and_analyze_specific_bacteria_gem{
	#输入是GEM文件，输出是所有的解析结果，包括代谢物，代谢反应，摄取物，分泌物
	#初始化，建立文件夹
	mkdir "tmp/BFS_result" unless -e "tmp/BFS_result"; 
	mkdir "tmp/BFS_result/storable" unless -e "tmp/BFS_result/storable";
	mkdir "tmp/BFS_result/metabolite" unless -e "tmp/BFS_result/metabolite";
	my %hash; #存储结果的哈希数据结构
	opendir Dir,$dir;
	my @files=readdir Dir; #读取目录
	@files=sort @files; 
	my %bac_names;
	#获取所有的代谢物
	foreach my $file(@files){ #遍历数组
		next if $file=~ /^\./; 
		$file=~ /^(.*?)\.xml/;  #检索每一个xml文件，用正则表达式获取GEM名称
		my $f_name=$1;
		next if $new_carve_me_not_bacteria{$f_name}; #bigg数据库里有一些人类、小鼠等非细菌的GEM
		$bac_names{$f_name}=1;
		print "$f_name\n";
		open IN,"<$dir/$file";
		open OUT,">tmp/BFS_result/metabolite/$f_name.txt";
		my $id;
		my $name;
		while(<IN>){ #逐行读取GEM文件
			last() if $_=~ /<\/listOfSpecies>/; #这里先识别代谢物的末XML节点
			next unless $_=~ /id="(.*?)".*?name="(.*?)"/ or $_=~ /id="(.*?)"/ or $_=~ /name="(.*?)"/; #正则表达式匹配，获取代谢物的所有信息
			if($_=~ / id="(.*?)".*?name="(.*?)"/){ #可能分为好几种情况
				$id=$1;
				$name=$2;
				print OUT "id=$id name=$name\n";
				$hash{$id}{"name"}=$name; #把信息存储到哈希数据结构中
			}elsif($_=~ /id="(.*?)"/){
				$id=$1;
			}elsif($_=~ /name="(.*?)"/){
				$name=$1;
				print OUT "id=$id name=$name\n";
				$hash{$id}{"name"}=$name;
			}
		}
		close IN;
		close OUT;
	}
	lock_nstore(\%bac_names,"tmp/BFS_result/storable/bacteria_names.storable"); #将数据结构存储为二进制格式文件
	
	open OUT,">tmp/BFS_result/all_metabolites.txt";
	foreach my $metabolite(keys %hash){
		print OUT "id=$metabolite name=$hash{$metabolite}{'name'}\n"; #将结果打印为txt文件
	}
	close OUT;
	lock_nstore(\%hash,"tmp/BFS_result/storable/metabolite.storable");
	open OUT,">tmp/BFS_result/all_metabolites_names.txt";
	my %names;
	foreach my $metabolite(keys %hash){
		$names{$hash{$metabolite}{'name'}}=$metabolite; #建立由代谢物ID到代谢物名称的哈希映射
	}
	foreach my $name(keys %names){
		print OUT "name=$name id=$names{$name}\n";
	}
	close OUT;
	lock_nstore(\%names,"tmp/BFS_result/storable/metabolite_name_to_id.storable");
	# return;

	#获取所有代谢反应
	my %all_reactions; #代谢反应哈希数据结构
	opendir Dir,$dir; #读取GEM目录
	my @strains=readdir Dir;
	@strains=sort @strains;
	my %metabolites=%{lock_retrieve("tmp/BFS_result/storable/metabolite.storable")}; #读取已经存储的代谢物哈希
	opendir Dir,$dir;
	@files=readdir Dir;
	@files=sort @files;
	$num=0;
	foreach my $file(@files){
		next if $file=~ /^\./;
		$file=~ /^(.*?)\.xml/;
		my $f_name=$1;
		next if $new_carve_me_not_bacteria{$f_name};
		print "analyzing $f_name\n";
		open IN,"<$dir/$file";
		my $isreaction=0; #这里的几个变量是为了判断XML里的代谢反应的反应物及产物节点
		my $isreactant=0;
		my $isproduct=0;
		my $reaction_id="";
		my $reaction_name="unknown"; #有些代谢反应没有名字，预设为unknown
		
		while(<IN>){ #逐行读取GEM文件
			last if $_=~ /<\/listOfReactions>/; #这里检索代谢反应的XML末节点，遇到这里就停止遍历
			if($_=~ /<reaction.*?id="(.*?)".*?name="(.*?)"/ or $_=~ /<reaction.*?id="(.*?)"/ or $_=~ /<reaction/ or $_=~ /<reaction.*?id="(.*?)".*?>/){
				if($_=~ /<reaction.*?id="(.*?)".*?name="(.*?)"/){ #用正则表达式匹配，也分为很多种情况
					$isreaction=1;
					$reaction_id=$1;
					$reaction_name=$2;
					$all_reactions{$reaction_id}{'reaction_name'}=$reaction_name;
				}elsif($_=~ /<reaction.*?id="(.*?)".*?>/){
					$isreaction=1;
					$reaction_id=$1;
					$all_reactions{$reaction_id}{'reaction_name'}=$reaction_name;
				}elsif($_=~ /<reaction.*?id="(.*?)"/){
					$isreaction=1;
					$reaction_id=$1;
				}elsif($_=~ /<reaction/){
					$isreaction=1;
				}
			}elsif($isreaction==1 and $_=~ /id="(.*?)".*?name="(.*?)"/){
				$reaction_id=$1;
				$reaction_name=$2;
				$all_reactions{$reaction_id}{'reaction_name'}=$reaction_name;
			}elsif($isreaction==1 and $_=~ /name="(.*?)"/){
				$reaction_name=$1;
				$all_reactions{$reaction_id}{'reaction_name'}=$reaction_name;
			}
			if($isreaction==1 and $_=~ /reversible="(.*?)"/){ #判断反应是否可逆
				my $reversible=$1;
				$all_reactions{$reaction_id}{'reversible'}=$reversible;
				if($reversible eq "true"){
					$num++; #这里记录可逆反应的数量
				}
			}
			$isreaction=0 if $_=~ /<\/reaction>/; #代谢反应节点已经结束
			$isreactant=1 if $_=~ /<listOfReactants>/; #反应物节点开始
			$isreactant=0 if $_=~ /<\/listOfReactants>/; #反应物节点结束
			$isproduct=1 if $_=~ /<listOfProducts>/; #产物节点开始
			$isproduct=0 if $_=~ /<\/listOfProducts>/; #产物节点结束
			if($_=~ /<speciesReference.*?species="(.*?)".*?stoichiometry="(.*?)".*?\/>/){ #用正则表达式匹配代谢反应里的反应物、产物信息
				my $species=$1;
				my $stoichiometry=$2;
				my $constant=$3;
				$all_reactions{$reaction_id}{'reactant'}{$species}{'id'}=$species if $isreactant==1; #如果是反应物，则加入对应哈希
				$all_reactions{$reaction_id}{'reactant'}{$species}{'name'}=$metabolites{$species}{'name'} if $isreactant==1;
				$all_reactions{$reaction_id}{'reactant'}{$species}{'stoichiometry'}=$stoichiometry if $isreactant==1;
				$all_reactions{$reaction_id}{'product'}{$species}{'id'}=$species if $isproduct==1; #如果是产物，则加入对应哈希
				$all_reactions{$reaction_id}{'product'}{$species}{'name'}=$metabolites{$species}{'name'} if $isproduct==1;
				$all_reactions{$reaction_id}{'product'}{$species}{'stoichiometry'}=$stoichiometry if $isproduct==1;
				if($all_reactions{$reaction_id}{'reversible'} eq "true"){ #如果反应可逆，则将该反应的逆向反应也加入进哈希中来
					my $rev_id="$reaction_id\_rev"; #命名为rev反应
					$all_reactions{$rev_id}{'reaction_name'}="$all_reactions{$reaction_id}{'reaction_name'}\_rev";
					$all_reactions{$rev_id}{'product'}{$species}{'id'}=$species if $isreactant==1;
					$all_reactions{$rev_id}{'product'}{$species}{'name'}=$metabolites{$species}{'name'} if $isreactant==1;
					$all_reactions{$rev_id}{'product'}{$species}{'stoichiometry'}=$stoichiometry if $isreactant==1;
					$all_reactions{$rev_id}{'reactant'}{$species}{'id'}=$species if $isproduct==1;
					$all_reactions{$rev_id}{'reactant'}{$species}{'name'}=$metabolites{$species}{'name'} if $isproduct==1;
					$all_reactions{$rev_id}{'reactant'}{$species}{'stoichiometry'}=$stoichiometry if $isproduct==1;
				}
			}
		}
		close IN;
	}
	#下面都是存储和打印哈希结构
	open OUT,">tmp/BFS_result/all_reactions.txt";
	print OUT Dumper(%all_reactions);
	close OUT;
	# return;
	open OUT,">tmp/BFS_result/all_reactions_statistics.txt";
	my $reaction_number=keys %all_reactions;
	print OUT "$reaction_number\n";
	print "reversible_reaction: $num\n";
	print "total_reaction: $reaction_number\n";
	# return;
	foreach my $reaction_id(sort keys %all_reactions){
		print OUT "reaction_id=$reaction_id ";
		print OUT "$all_reactions{$reaction_id}{'reaction_name'}\n";
		print OUT "reactant=";
		foreach my $reactant(sort keys %{$all_reactions{$reaction_id}{'reactant'}}){
			print OUT "$all_reactions{$reaction_id}{'reactant'}{$reactant}{'name'} ";
		}
		print OUT "\n";
		print OUT "product=";
		foreach my $product(sort keys %{$all_reactions{$reaction_id}{'product'}}){
			print OUT "$all_reactions{$reaction_id}{'product'}{$product}{'name'} ";
		}
		print OUT "\n";
	}
	close OUT;
	
	lock_nstore(\%all_reactions,"tmp/BFS_result/storable/all_reactions.storable");
	# return;

	#获取每一个细菌GEM的代谢反应，并分细菌存储
	mkdir "tmp/BFS_result/reaction_of_every_bacteria" unless -e "tmp/BFS_result/reaction_of_every_bacteria";
	mkdir "tmp/BFS_result/reaction_of_every_bacteria/storable" unless -e "tmp/BFS_result/reaction_of_every_bacteria/storable";
	my %bac_reactions; #每一个细菌的代谢反应哈希
	opendir Dir,$dir;
	@files=readdir Dir;
	@files=sort @files;
	foreach my $file(@files){
		next if $file=~ /^\./;
		$file=~ /^(.*?)\.xml/;
		my $f_name=$1;
		next if $new_carve_me_not_bacteria{$f_name};
		open IN,"<$dir/$file";
		my $reaction_id="";
		my $isreaction=0;
		while(<IN>){ #逐行读取每一个细菌的GEM文件
			last if $_=~ /<\/listOfReactions>/;
			if($_=~ /<reaction.*?id="(.*?)"/ or $_=~ /<reaction/){
				if($_=~ /<reaction.*?id="(.*?)"/){
					$isreaction=1;
					$reaction_id=$1;
					$bac_reactions{$f_name}{$reaction_id}=1;
					if($all_reactions{$reaction_id}{'reversible'} eq "true"){
						my $rev_id="$reaction_id\_rev";
						$bac_reactions{$f_name}{$rev_id}=1;
					}
				}elsif($_=~ /<reaction/){
					$isreaction=1;
				}
			}elsif($isreaction==1 and $_=~ /id="(.*?)"/){
				$reaction_id=$1;
				$bac_reactions{$f_name}{$reaction_id}=1;
				if($all_reactions{$reaction_id}{'reversible'} eq "true"){
					my $rev_id="$reaction_id\_rev";
					$bac_reactions{$f_name}{$rev_id}=1;
				}
			}
		}
		close IN;
		my $reaction_number=keys %{$bac_reactions{$f_name}};
		print "reaction_number=$reaction_number\n";
		open OUT,">tmp/BFS_result/reaction_of_every_bacteria/$f_name.txt";
		print OUT Dumper(%{$bac_reactions{$f_name}});
		close OUT;
		lock_nstore(\%{$bac_reactions{$f_name}},"tmp/BFS_result/reaction_of_every_bacteria/storable/$f_name.storable");
	}
	# return;
	foreach my $bacteria(sort keys %bac_names){
		print "get reaction info for each bacteria: $bacteria\n";
		mkdir "tmp/BFS_result/reaction_info_for_every_bacteria" unless -e "tmp/BFS_result/reaction_info_for_every_bacteria";
		mkdir "tmp/BFS_result/reaction_info_for_every_bacteria/storable" unless -e "tmp/BFS_result/reaction_info_for_every_bacteria/storable";
		mkdir "tmp/BFS_result/reaction_info_for_every_bacteria/reaction_info" unless -e "tmp/BFS_result/reaction_info_for_every_bacteria/reaction_info";
		my %bac_reaction_id=%{lock_retrieve("tmp/BFS_result/reaction_of_every_bacteria/storable/$bacteria.storable")};
		my %bac;
		foreach my $r_id(keys %bac_reaction_id){
			%{$bac{$r_id}}=%{$all_reactions{$r_id}};
		}
		open OUT,">tmp/BFS_result/reaction_info_for_every_bacteria/reaction_info/$bacteria.txt";
		print OUT Dumper(%bac);
		close OUT;
		lock_nstore(\%bac,"tmp/BFS_result/reaction_info_for_every_bacteria/storable/$bacteria.storable");
	}
	# return;

	#获取每一个代谢反应的产物、反应物，各仅存储一个相关反应 
	#get_reaction
	my %all_reactants;
	my %all_products;
	mkdir "tmp/BFS_result/reaction" unless -e "tmp/BFS_result/reaction";
	mkdir "tmp/BFS_result/reaction/reactant" unless -e "tmp/BFS_result/reaction/reactant";
	mkdir "tmp/BFS_result/reaction/product" unless -e "tmp/BFS_result/reaction/product";
	opendir Dir,$dir;
	@files=readdir Dir;
	@files=sort @files;
	foreach my $file(@files){
		next if $file=~ /^\./;
		$file=~ /^(.*?)\.xml/;
		my $f_name=$1;
		next if $new_carve_me_not_bacteria{$f_name};
		print "get reactant and product for each bacteria: $f_name\n";
		my %bac_reaction=%{lock_retrieve("tmp/BFS_result/reaction_info_for_every_bacteria/storable/$f_name.storable")};
		my %reactants;
		my %products;
		open OUT1,">tmp/BFS_result/reaction/reactant/$f_name.txt";
		open OUT2,">tmp/BFS_result/reaction/product/$f_name.txt";
		foreach my $reaction_id(sort keys %bac_reaction){
			foreach my $reactant(sort keys %{$bac_reaction{$reaction_id}{'reactant'}}){
				# 判断产物和摄取物时需要排除EX反应！
				next if $reaction_id=~ /R_EX_/;
				$reactants{$reactant}{'name'}=$bac_reaction{$reaction_id}{'reactant'}{$reactant}{'name'};
				$reactants{$reactant}{'reaction_id'}=$reaction_id;
				$all_reactants{$reactant}{'name'}=$bac_reaction{$reaction_id}{'reactant'}{$reactant}{'name'};
				$all_reactants{$reactant}{'reaction_id'}=$reaction_id;
			}
			foreach my $product(sort keys %{$bac_reaction{$reaction_id}{'product'}}){
				# 判断产物和摄取物时需要排除EX反应！
				next if $reaction_id=~ /R_EX_/;
				$products{$product}{'name'}=$bac_reaction{$reaction_id}{'product'}{$product}{'name'};
				$products{$product}{'reaction_id'}=$reaction_id;
				$all_products{$product}{'name'}=$bac_reaction{$reaction_id}{'product'}{$product}{'name'};
				$all_products{$product}{'reaction_id'}=$reaction_id;
			}
		}
		foreach my $reactant(keys %reactants){
			print OUT1 "id=$reactant name=$reactants{$reactant}{'name'} reaction_id=$reactants{$reactant}{'reaction_id'}\n";
		}
		foreach my $product(keys %products){
			print OUT2 "id=$product name=$products{$product}{'name'} reaction_id=$products{$product}{'reaction_id'}\n";
		}
		close OUT1;
		close OUT2;
	}
	open OUT1,">tmp/BFS_result/all_reactants.txt";
	open OUT2,">tmp/BFS_result/all_products.txt";
	foreach my $reactant(keys %all_reactants){
		print OUT1 "id=$reactant name=$all_reactants{$reactant}{'name'} reaction_id=$all_reactants{$reactant}{'reaction_id'}\n";
	}
	foreach my $product(keys %all_products){
		print OUT2 "id=$product name=$all_products{$product}{'name'} reaction_id=$all_products{$product}{'reaction_id'}\n";
	}
	lock_nstore(\%all_reactants,"tmp/BFS_result/storable/reactants.storable");
	lock_nstore(\%all_products,"tmp/BFS_result/storable/products.storable");
	# return;

	#获取所有的摄取物
	#get_intake
	mkdir "tmp/BFS_result/intake" unless -e "tmp/BFS_result/intake";
	mkdir "tmp/BFS_result/secretion" unless -e "tmp/BFS_result/secretion";
	my $dir2="tmp/BFS_result/reaction/reactant";
	opendir DIR2,$dir2;
	@files=readdir DIR2;
	my %all_intakes;
	# return;
	foreach my $file(@files){
		next if $file=~ /^\./;
		$file=~ /^(.*?)\.txt/;
		my $f_name=$1;
		next if $new_carve_me_not_bacteria{$f_name};
		next unless ($bac_names{$f_name} and $bac_names{$f_name}==1);
		my %intakes;
		open IN,"<tmp/BFS_result/reaction/reactant/$file"; #遍历已经找到的所有代谢反应的反应物信息
		while(<IN>){
			chomp;
			$_=~ /id=(.*?) name=(.*?) reaction_id=(.*?)$/;
			my $id=$1;
			my $name=$2;
			print "$id $name\n";
			$intakes{$metabolites{$id}{'name'}}="" if $id=~ /_e$/; #如果反应物是在胞外的，则说明其是一种会被摄取到胞内的物质
			$all_intakes{$metabolites{$id}{'name'}}{'bacteria'}{$f_name}="" if $id=~ /_e$/;
		}
		$file=~ /^(.*?)\.txt/;
		open OUT,">tmp/BFS_result/intake/$1.txt";
		foreach my $intake(keys %intakes){
			print OUT "$intake\n";
		}
		close OUT;
	}
	open OUT,">tmp/BFS_result/all_intake_names.txt";
	foreach my $intake(keys %all_intakes){
		print OUT "$intake\n";
	}
	close OUT;
	open OUT,">tmp/BFS_result/all_intake_with_bacteria.txt";
	# print Dumper(%all_intakes);
	close OUT;
	lock_nstore(\%all_intakes,"tmp/BFS_result/storable/all_intake_with_bacteria.storable");
	# return;

	#获取所有的分泌物
	#get_secretion
	my $dir3="tmp/BFS_result/reaction/product";
	opendir DIR3,$dir3;
	@files=readdir DIR3;
	my %all_secretions;
	foreach my $file(@files){
		next if $file=~ /^\./;
		$file=~ /^(.*?)\.txt/;
		my $f_name=$1;
		next if $new_carve_me_not_bacteria{$f_name};
		next unless ($bac_names{$f_name} and $bac_names{$f_name}==1);
		print "$file\n";
		my %bac_reaction=%{lock_retrieve("tmp/BFS_result/reaction_info_for_every_bacteria/storable/$f_name.storable")};
		my %secretions;
		my %real_secretions;
		open IN,"<tmp/BFS_result/reaction/product/$file";
		while(<IN>){ #遍历的是已经找到的所有代谢反应的产物
			chomp;
			$_=~ /id=(.*?) name=(.*?) reaction_id=(.*?)$/;
			my $id=$1;
			my $name=$2;
			my $reaction_id=$3;
			$secretions{$id}{'name'}=$name;
			$secretions{$id}{'reaction_id'}=$reaction_id;
		}
		foreach my $id(keys %secretions){
			my $name=$secretions{$id}{'name'};
			my $reaction_id=$secretions{$id}{'reaction_id'};
			next unless $id=~ /_e$/; #只有该产物是胞外的，其才可能是一种分泌物，但是还分为以下几种情况
			foreach my $reactant(keys %{$bac_reaction{$reaction_id}{'reactant'}}){ #要分别看其对应的代谢反应的产物是在哪个亚细胞定位
				next unless $bac_reaction{$reaction_id}{'reactant'}{$reactant}{'name'};
				my $r_name=$bac_reaction{$reaction_id}{'reactant'}{$reactant}{'name'};
				next unless $r_name eq $name;
				if($reactant=~ /_c$/){ #如果是从c直接到e，则其一定是分泌物
					$real_secretions{$name}="";
					$all_secretions{$name}{'bacteria'}{$f_name}=""
				}elsif($reactant=~ /_p$/){ #如果是从p到e，则要再看一下有没有从c到p的反应
					my $n_reaction_id=$secretions{$reactant}{'reaction_id'}; #继续找到p物质作为产物对应的代谢反应
					next unless $bac_reaction{$n_reaction_id}{'reactant'};
					foreach my $n_reactant(keys %{$bac_reaction{$n_reaction_id}{'reactant'}}){
						next unless $bac_reaction{$n_reaction_id}{'reactant'}{$n_reactant}{'name'};
						my $n_r_name=$bac_reaction{$n_reaction_id}{'reactant'}{$n_reactant}{'name'};
						
						next unless $n_r_name eq $name;
						if($n_reactant=~ /_c$/){ #如果p物质作为产物参与的代谢反应的反应物为c，则证明存在c到p的反应，则其为分泌物
							$real_secretions{$name}="";
							$all_secretions{$name}{'bacteria'}{$f_name}="";
						}
					}
				}
			}
		}
		$file=~ /^(.*?)\.txt/;
		open OUT,">tmp/BFS_result/secretion/$1.txt";
		foreach my $secretion(keys %real_secretions){
			print OUT "$secretion\n";
		}
		close OUT;
		my $secretion_number=keys %real_secretions; #获取细菌的所有分泌物
		print "$1 secretion_number=$secretion_number\n";
	}
	# return;
	open OUT,">tmp/BFS_result/all_secretion_names.txt";
	foreach my $secretion(keys %all_secretions){
		print OUT "$secretion\n";
	}
	close OUT;
	open OUT,">tmp/BFS_result/all_secretion_with_bacteria.txt";
	print OUT Dumper(%all_secretions);
	close OUT;
	lock_nstore(\%all_secretions,"tmp/BFS_result/storable/all_secretion_with_bacteria.storable");
	# return;

	#获取每一个细菌的所有摄取物和分泌物
	#get_intake_and_secretion_of_every_bacteria
	mkdir "tmp/BFS_result/intake_and_secretion_of_every_bacteria" unless -e "tmp/BFS_result/intake_and_secretion_of_every_bacteria";
	mkdir "tmp/BFS_result/intake_and_secretion_of_every_bacteria/storable" unless -e "tmp/BFS_result/intake_and_secretion_of_every_bacteria/storable";
	my %bac;
	my $product_dir="tmp/BFS_result/reaction/product";
	opendir DIR,$product_dir;
	@files=readdir DIR;
	my $num=0;
	foreach my $file(sort @files){
		next if $file=~ /^\./;
		$num++;
		my %secretions;
		open IN,"<tmp/BFS_result/reaction/product/$file";
		$file=~ /^(.*?)\.txt/;
		my $bac_name=$1;
		next if $new_carve_me_not_bacteria{$bac_name};
		while(<IN>){
			chomp;
			$_=~ /id=(.*?) name=(.*?) reaction_id=(.*?)$/;
			my $id=$1;
			my $name=$2;
			my $reaction_id=$3;
			# 下面有bug，把M_e开头的代谢物也计入了
			$bac{$bac_name}{'secretion'}{$id}{'name'}=$name if $id=~ /_e$/; #把之前已经获取的信息按照细菌分别存储起来
			$bac{$bac_name}{'secretion'}{$id}{'reaction_id'}=$reaction_id if $id=~ /_e$/;
		}
	}
	closedir DIR;
	$num=0;
	my $reactant_dir="tmp/BFS_result/reaction/reactant";
	opendir DIR,$reactant_dir;
	@files=readdir DIR;
	foreach my $file(sort @files){
		next if $file=~ /^\./;
		$num++;
		print "$num\n" if $num%100==0;
		my %intakes;
		open IN,"<tmp/BFS_result/reaction/reactant/$file";
		$file=~ /^(.*?)\.txt/;
		my $bac_name=$1;
		next if $new_carve_me_not_bacteria{$bac_name};
		while(<IN>){
			chomp;
			$_=~ /id=(.*?) name=(.*?) reaction_id=(.*?)$/;
			my $id=$1;
			my $name=$2;
			my $reaction_id=$3;
			$bac{$bac_name}{'intake'}{$id}{'name'}=$name if $id=~ /_e$/;
			$bac{$bac_name}{'intake'}{$id}{'reaction_id'}=$reaction_id if $id=~ /_e$/;
		}
	}
	open OUT,">tmp/BFS_result/intake_and_secretion_of_every_bacteria.txt";
	print OUT Dumper(%bac);
	close OUT;
	lock_nstore(\%bac,"tmp/BFS_result/storable/intake_and_secretion_of_every_bacteria.storable");
	foreach my $bacteria(sort keys %bac){
		open OUT,">tmp/BFS_result/intake_and_secretion_of_every_bacteria/$bacteria.txt";
		print OUT Dumper(%{$bac{$bacteria}});
		close OUT;
		lock_nstore(\%{$bac{$bacteria}},"tmp/BFS_result/intake_and_secretion_of_every_bacteria/storable/$bacteria.storable");
	}
	# return;
}


#寻找从摄取物到分泌物的代谢通路
sub new_carve_me_find_path_from_secretion_to_intake_for_each_bacteria{
	my %bac_names=%{lock_retrieve("tmp/BFS_result/storable/bacteria_names.storable")};
	mkdir "tmp/BFS_result/path_from_secretion_to_intake" unless -e "tmp/BFS_result/path_from_secretion_to_intake";
	my %path;
	my %path_intake_to_media;

	#初始化，构建代谢反应的产物-反应物的双向映射
	foreach my $bac(sort keys %bac_names){
		next if $bac eq "Recon3D";
		print "$bac\n";
		my %bac_reaction=%{lock_retrieve("tmp/BFS_result/reaction_info_for_every_bacteria/storable/$bac.storable")}; #读取细菌代谢反应哈希
		my %product_to_reactant; 
		my %bac_metabolite=%{lock_retrieve("tmp/BFS_result/intake_and_secretion_of_every_bacteria/storable/$bac.storable")}; #读取细菌代谢物哈希
		# return;
		my %bac_meta_name; #获取细菌的所有代谢物名称和id
		foreach my $meta_id(keys %bac_metabolite){
			next unless $bac_metabolite{$meta_id}{'name'};
			my $meta_name=$bac_metabolite{$meta_id}{'name'};
			$bac_meta_name{$meta_name}=1;
		}
		my %bac_intake=%{$bac_metabolite{'intake'}}; #获取细菌的所有摄取物
		my %bac_int_name;
		foreach my $intake(keys %bac_intake){
			my $intake_name=$bac_intake{$intake}{'name'};
			$bac_int_name{$intake_name}=1;
		} 
		my %bac_secretion=%{$bac_metabolite{'secretion'}};
		# return;
		my %sec_name;
		my $secretion_number=keys %bac_secretion; #获取细菌的所有分泌物
		print "$bac secretion_number=$secretion_number\n";
		foreach my $secretion(keys %bac_secretion){
			my $secretion_name=$bac_secretion{$secretion}{'name'};
			$sec_name{$secretion_name}=1;
		} 
		open OUT,">tmp/BFS_result/path_from_secretion_to_intake/$bac.txt";
		foreach my $reaction_id(keys %bac_reaction){ #建立由代谢反应产物->代谢反应反应物的映射
			foreach my $product(keys %{$bac_reaction{$reaction_id}{'product'}}){
				next unless $bac_reaction{$reaction_id}{'product'}{$product}{'name'};
				my $product_name=$bac_reaction{$reaction_id}{'product'}{$product}{'name'};
				next if $basic_molecules{$product_name};
				foreach my $reactant(keys %{$bac_reaction{$reaction_id}{'reactant'}}){
					my $reactant_name=$bac_reaction{$reaction_id}{'reactant'}{$reactant}{'name'};
					next if $basic_molecules{$reactant_name};
					next if $product_name eq $reactant_name;
					$product_to_reactant{$product_name}{'reactant'}{$reactant_name}{'reaction'}{$reaction_id}="";
				}
			}
		}
		foreach my $product(keys %add_reaction){
			$product_to_reactant{$product}{'reactant'}{$add_reaction{$product}}{'reaction'}{'add_reaction'}="";
		}
		# return;
		
		#利用分层遍历，寻找从摄取物到培养基底物的代谢通路
		#目的是判断一些摄取物是否可以用培养基底物合成
		foreach my $intake(sort keys %bac_int_name){
			next if $basic_molecules{$intake};
			my @arr;
			push @arr,$intake;
			my %circle;
			my $level=-1;
			while(@arr){
				my $size=$#arr+1;
				last unless $arr[0];
				$level++;
				last if $level==$max_level; #####################################
				while($size--){
					my $t=shift @arr;
					next if $circle{$t} and $circle{$t}==1;
					$circle{$t}=1;
					if($t ne $intake and $target_media{$t} and $target_media{$t}==1){
						$path_intake_to_media{$bac}{$intake}{'media'}{$t}{'level'}=$level;
						next;
					}
					next unless $product_to_reactant{$t};
					foreach my $reactant(keys %{$product_to_reactant{$t}{'reactant'}}){
						push @arr,$reactant;
					}
				}
				@arr=unique(@arr);
			}
		}
		# return;

		#利用分层遍历，寻找从分泌物到摄取物的代谢通路
		foreach my $secretion(sort keys %bac_secretion){ #每一轮遍历的起点是细菌的分泌物
			$secretion=$bac_secretion{$secretion}{'name'}; 
			next if $basic_molecules{$secretion}; #基础小分子物质不参与遍历
			my @arr;  #用于存储遍历队列的数组
			push @arr,$secretion; #遍历起点
			my %circle; #用于防治环形遍历路径的辅助变量
			my $level=-1; #分层遍历的层数
			while(@arr){ #当遍历队列非空
				my $size=$#arr+1; 
				last unless $arr[0];
				$level++;
				last if $level==$max_level; #############################################
				while($size--){
					my $t=shift @arr;
					next if $circle{$t} and $circle{$t}==1; #防止环形路径产生
					$circle{$t}=1; 
					if($t ne $secretion and $bac_int_name{$t} and $bac_int_name{$t}==1){ #如果遍历到了摄取物，而且摄取物与分泌物不同
						$path{$bac}{$secretion}{'intake'}{$t}{'level'}=$level; #记录路径
						if($path_intake_to_media{$bac}{$t}){ 
							foreach my $media(sort keys %{$path_intake_to_media{$bac}{$t}{'media'}}){
								$path{$bac}{$secretion}{'intake'}{$t}{'media'}{$media}=$path_intake_to_media{$bac}{$t}{'media'}{$media}{'level'};
							}
						}
						print OUT "$secretion linked to $t with level=$level\n";
						next;
					}
					next unless $product_to_reactant{$t}; #如果已经没有以该物质作为产物的代谢反应，则进入下一轮循环；否则，将其作为产物对应的代谢反应的反应物加入到遍历队列中
					foreach my $reactant(keys %{$product_to_reactant{$t}{'reactant'}}){
						push @arr,$reactant; 
					}
				}
				@arr=unique(@arr); #遍历队列去重
			}
			print OUT "\n";
		}
		close OUT;
	}
	open OUT,">tmp/BFS_result/path_from_secretion_to_intake.txt";
	print OUT Dumper(%path);
	close OUT;
	open OUT,">tmp/BFS_result/path_from_intake_to_media.txt";
	print OUT Dumper(%path_intake_to_media);
	close OUT;
	lock_nstore(\%path,"tmp/BFS_result/storable/path_from_secretion_to_intake.storable");
	lock_nstore(\%path_intake_to_media,"tmp/BFS_result/storable/path_from_intake_to_media.storable");
}

#按照代谢目标，筛选细菌组合
sub new_carve_me_analyze_path_result{
	my $target=shift @_; #目标物质
	my $intermediate=shift @_; #代谢衔接物
	my $csource=shift @_; #碳源
	my %path=%{lock_retrieve("tmp/BFS_result/storable/path_from_secretion_to_intake.storable")}; #之前对每一种细菌获取的分泌物到摄取物的路径
	my %result;
	#my %bac_id_to_name=%{lock_retrieve("tmp/BFS_result/storable/bac_id_to_name.storable")}; 
	$result{'target'}{'name'}=$target;
	$result{'intermediate'}{'name'}=$intermediate;
	#寻找从代谢衔接物到目标物质的细菌
	foreach my $bac(sort keys %path){ #路径内到每一种细菌
		my $bac_name=$bac;
		next unless $path{$bac}{$target};
		foreach my $intake(sort keys %{$path{$bac}{$target}{'intake'}}){ #该细菌该目标物质对应的摄取物
			if($intake eq $intermediate){
				$result{'target'}{'bacteria'}{$bac_name}=$bac;
			}
			if($path{$bac}{$target}{'intake'}{$intake}{'media'}){ #判断该摄取物是否可以由对应的培养基底物直接合成
				foreach my $media(keys %{$path{$bac}{$target}{'intake'}{$intake}{'media'}}){
					if($media eq $intermediate){
						$result{'target'}{'bacteria'}{$bac_name}=$bac;
					}
				}
			}
		}
	}
	#寻找从底物到代谢衔接物的细菌
	foreach my $bac(sort keys %path){
		my $bac_name=$bac;
		next unless $path{$bac}{$intermediate};
		foreach my $intake(sort keys %{$path{$bac}{$intermediate}{'intake'}}){
			if($intake eq $csource){
				$result{'intermediate'}{'bacteria'}{$bac_name}=$bac;
			}
			if($path{$bac}{$intermediate}{'intake'}{$intake}{'media'}){
				foreach my $media(keys %{$path{$bac}{$intermediate}{'intake'}{$intake}{'media'}}){
					if($media eq $csource){
						$result{'intermediate'}{'bacteria'}{$bac_name}=$bac;
					}
				}
			}
		}
	}
	print Dumper(%result);
	open OUT,">tmp/BFS_result/$target\_$intermediate\_$csource.txt";
	print OUT Dumper(%result);
}

sub new_carve_me_find_all_secretion_can_be_synthesized_from_specific_intermediate_csource{ #寻找可以上述物质为碳源合成的所有产物及相应细菌
	my @csources=@{shift @_}; #获取碳源列表
	my %path;
	%path=%{lock_retrieve("tmp/BFS_result/storable/path_from_secretion_to_intake.storable")}; #从分泌物到摄取物的路径及可以实现转化的细菌列表
	mkdir "tmp/BFS_result/all_secretion_can_be_synthesized_from_specific_intermediate_csource" unless -e "tmp/BFS_result/all_secretion_can_be_synthesized_from_specific_intermediate_csource";
	my %intake_with_bac=%{lock_retrieve("tmp/BFS_result/storable/all_intake_with_bacteria.storable")}; #可以摄取特定物质的细菌列表
	foreach my $csource(@csources){ #对每一个碳源分别遍历
		my %result;
		print "$csource\n";
		open OUT,">tmp/BFS_result/all_secretion_can_be_synthesized_from_specific_intermediate_csource/$csource.txt";
		open OUT2,">tmp/BFS_result/all_secretion_can_be_synthesized_from_specific_intermediate_csource/$csource\_with_bacteria.txt";
		foreach my $bac(sort keys %{$intake_with_bac{$csource}{'bacteria'}}){ #遍历所有可以摄取该物质的细菌
			foreach my $target(keys %{$path{$bac}}){ #该细菌可以合成的所有目标物质
				foreach my $intake(sort keys %{$path{$bac}{$target}{'intake'}}){ #如果该细菌合成该目标物质所使用的碳源与列表中的碳源一致
					if($intake eq $csource){
						$result{$target}{'bacteria'}{$bac}=1;
					}
				}
			}
		}
		# print Dumper(%result);
		print OUT2 Dumper(%result);
		foreach my $target(sort keys %result){
			print OUT "$target\n";
		}
		close OUT;
		close OUT2;
		open OUT,">tmp/BFS_result/all_secretion_can_be_synthesized_from_specific_intermediate_csource/$csource.csv";
		foreach my $target(sort keys %result){
			print OUT "$target\n";
			foreach my $bac(sort keys %{$result{$target}{'bacteria'}}){
				print OUT "$bac\n";
			}
			print OUT "\n\n";
		}
		close OUT;
	}
}

sub new_carve_me_find_bacteria_can_produce_with_csource{ #寻找可以用特定碳源合成上述物质的细菌
	my @targets=@{shift @_}; #目标物质
	my $csource=shift @_; #碳源
	my %path=%{lock_retrieve("tmp/BFS_result/storable/path_from_secretion_to_intake.storable")}; #从分泌物到摄取物的路径及可以实现转化的细菌列表
	my %result;
	my %bac_id_to_name=%{lock_retrieve("tmp/BFS_result/storable/bac_id_to_name.storable")};
	foreach my $target(@targets){ #遍历目标物质
		foreach my $bac(sort keys %path){ #遍历路径中的每一种细菌
			my $bac_name=$bac;
			next unless $path{$bac}{$target};
			foreach my $intake(sort keys %{$path{$bac}{$target}{'intake'}}){
				if($intake eq $csource){ #如果该细菌可以利用碳源合成该目标物质
					$result{$target}{'bacteria'}{$bac_name}=$bac;
				}
				if($path{$bac}{$target}{'intake'}{$intake}{'media'}){
					foreach my $media(keys %{$path{$bac}{$target}{'intake'}{$intake}{'media'}}){
						if($media eq $csource){
							$result{'target'}{'bacteria'}{$bac_name}=$bac;
						}
					}
				}
			}
		}
	}
	print Dumper(%result);
}

sub new_carve_me_deep_analysis_of_bacteria_internal_pathways_and_link_target_secretion_to_target_intake{ #详细的查看对应细菌内由某物质合成某物质涉及的代谢通路
	my $target_intake=shift @_; #获取目标碳源
	my $target_secretion=shift @_; #获取目标物质
	my $bac_name=shift @_; #获取细菌名称
	print "$target_intake\n";
	print "$target_secretion\n";
	print "$bac_name\n";
	my %metabolite=%{lock_retrieve("tmp/BFS_result/storable/metabolite.storable")};
	my %secretion_and_intake=%{lock_retrieve("tmp/BFS_result/intake_and_secretion_of_every_bacteria/storable/$bac_name.storable")};
	my %intake=%{$secretion_and_intake{'intake'}};
	foreach my $intake(keys %intake){
		$intake{$metabolite{$intake}{'name'}}=1;
	}
	mkdir "tmp/BFS_result/reactant_and_product" unless -e "tmp/BFS_result/reactant_and_product";
	mkdir "tmp/BFS_result/reactant_and_product/storable" unless -e "tmp/BFS_result/reactant_and_product/storable";
	my %bac_reaction=%{lock_retrieve("tmp/BFS_result/reaction_info_for_every_bacteria/storable/$bac_name.storable")};
	my %reactant_to_product;
	foreach my $reaction_id(keys %bac_reaction){
		my $reaction_name=$bac_reaction{$reaction_id}{'reaction_name'};
		foreach my $reactant_id(keys %{$bac_reaction{$reaction_id}{'reactant'}}){
			my $reactant_number=keys %{$bac_reaction{$reaction_id}{'reactant'}};
			my $product_number=keys %{$bac_reaction{$reaction_id}{'product'}};
			next if $product_number==0;
			if($reactant_number==1 and $product_number==1){
				my $r_name="";
				my $p_name="";
				foreach my $reactant(sort keys %{$bac_reaction{$reaction_id}{'reactant'}}){
					$r_name=$bac_reaction{$reaction_id}{'reactant'}{$reactant}{'name'};
				}
				foreach my $product(sort keys %{$bac_reaction{$reaction_id}{'product'}}){
					$p_name=$bac_reaction{$reaction_id}{'product'}{$product}{'name'};
				}
				next if $r_name eq $p_name;
			}
			if($reactant_number==$product_number){
				my @r_arr;
				my @p_arr;
				foreach my $reactant(sort keys %{$bac_reaction{$reaction_id}{'reactant'}}){
					push @r_arr,$bac_reaction{$reaction_id}{'reactant'}{$reactant}{'name'};
				}
				foreach my $product(sort keys %{$bac_reaction{$reaction_id}{'product'}}){
					push @p_arr,$bac_reaction{$reaction_id}{'product'}{$product}{'name'};
				}
				my $same=same(\@r_arr,\@p_arr);
				next if $same==1;
			}
			my $reactant_name=$bac_reaction{$reaction_id}{'reactant'}{$reactant_id}{'name'};
			foreach my $product(sort keys %{$bac_reaction{$reaction_id}{'product'}}){
				my $p_name=$bac_reaction{$reaction_id}{'product'}{$product}{'name'};
				next if $basic_molecules{$p_name};
				$reactant_to_product{$reactant_name}{'product'}{$p_name}{'reaction'}{$reaction_id}="";
			}
		}
	}
	foreach my $product(keys %add_reaction){
		$reactant_to_product{$add_reaction{$product}}{'product'}{$product}{'reaction'}{'add_reaction'}="";
	}
	open OUT,">tmp/BFS_result/reactant_and_product/$bac_name\_link_reactant_to_product.txt";
	print OUT Dumper(%reactant_to_product);
	close OUT;
	lock_nstore(\%reactant_to_product,"tmp/BFS_result/reactant_and_product/storable/$bac_name\_link_reactant_to_product.storable");
	my %product_to_reactant;
	foreach my $reaction_id(sort keys %bac_reaction){
		my $reactant_number=keys %{$bac_reaction{$reaction_id}{'reactant'}};
		my $product_number=keys %{$bac_reaction{$reaction_id}{'product'}};
		next if $product_number==0;
		if($reactant_number==1 and $product_number==1){
			my $r_name="";
			my $p_name="";
			foreach my $reactant(sort keys %{$bac_reaction{$reaction_id}{'reactant'}}){
				$r_name=$bac_reaction{$reaction_id}{'reactant'}{$reactant}{'name'};
			}
			foreach my $product(sort keys %{$bac_reaction{$reaction_id}{'product'}}){
				$p_name=$bac_reaction{$reaction_id}{'product'}{$product}{'name'};
			}
			next if $r_name eq $p_name;
		}
		if($reactant_number==$product_number){
			my @r_arr;
			my @p_arr;
			foreach my $reactant(sort keys %{$bac_reaction{$reaction_id}{'reactant'}}){
				push @r_arr,$bac_reaction{$reaction_id}{'reactant'}{$reactant}{'name'};
			}
			foreach my $product(sort keys %{$bac_reaction{$reaction_id}{'product'}}){
				push @p_arr,$bac_reaction{$reaction_id}{'product'}{$product}{'name'};
			}
			my $same=same(\@r_arr,\@p_arr);
			next if $same==1;
		}
		foreach my $product(sort keys %{$bac_reaction{$reaction_id}{'product'}}){
			my $p_name=$bac_reaction{$reaction_id}{'product'}{$product}{'name'};
			next if $p_name eq "H+" or $p_name eq "H2O" or $p_name eq "ADP";
			foreach my $reactant(sort keys %{$bac_reaction{$reaction_id}{'reactant'}}){
				my $r_name=$bac_reaction{$reaction_id}{'reactant'}{$reactant}{'name'};
				next if $basic_molecules{$r_name};
				$product_to_reactant{$p_name}{'reactant'}{$r_name}{'reaction'}{$reaction_id}="";
			}
		}
	}
	foreach my $product(keys %add_reaction){
		$product_to_reactant{$product}{'reactant'}{$add_reaction{$product}}{'reaction'}{'add_reaction'}="";
	}
	open OUT,">tmp/BFS_result/reactant_and_product/$bac_name\_link_product_to_reactant.txt";
	print OUT Dumper(%product_to_reactant);
	close OUT;
	lock_nstore(\%product_to_reactant,"tmp/BFS_result/reactant_and_product/storable/$bac_name\_link_product_to_reactant.storable");

	my @arr;
	my %path;
	push @{$path{$target_intake}{'metabolite'}},"$target_intake";
	push @{$path{$target_intake}{'reaction'}},"$target_intake";
	push @arr,$target_intake;
	my $level=-1;
	my $stop=0;
	my %circle;
	while(@arr){
		my $size=$#arr+1;
		$level++;
		print "\nlevel: $level\n";
		last if $stop==1;
		last if $level==$max_level; #################################
		while($size--){
			my $t=shift @arr;
			next if $circle{$t} and $circle{$t}==1;
			$circle{$t}=1;
			next if $basic_molecules{$t} and $t ne $target_intake;
			print "$t\n";
			foreach my $product(keys %{$reactant_to_product{$t}{'product'}}){
				foreach my $metabolite_path(sort @{$path{$t}{'metabolite'}}){
					push @{$path{$product}{'metabolite'}},"$metabolite_path->$product";
				}
				@{$path{$product}{'metabolite'}}=unique(@{$path{$product}{'metabolite'}});
				foreach my $reaction_path(sort @{$path{$t}{'reaction'}}){
					foreach my $reaction(keys %{$reactant_to_product{$t}{'product'}{$product}{'reaction'}}){
						push @{$path{$product}{'reaction'}},"$reaction_path->$reaction->$product";
					}
				}
				
				push @arr,$product;
				$stop=1 if $product eq $target_secretion;
			}
		}
		@{$path{$target_secretion}{'reaction'}}=unique(@{$path{$target_secretion}{'reaction'}}) if $path{$target_secretion}{'reaction'};
		@arr=unique(@arr);
	}
	print Dumper(@{$path{$target_secretion}{'metabolite'}});
	print Dumper(@{$path{$target_secretion}{'reaction'}});
	print "\n\n";
	return;
	foreach my $reaction_path(@{$path{$target_secretion}{'reaction'}}){
		print "$reaction_path:\n";
		my @reactions=split(/->/,$reaction_path);
		my $num=1;
		foreach my $reaction(@reactions){
			next unless $reaction=~ /^R_/;
			print "$num: ";
			$num++;
			print_reaction($reaction,\%bac_reaction);
		}
		print "\n\n";
	}
} 

sub same{
	my @arr1=@{shift @_};
	my @arr2=@{shift @_};
	@arr1=sort @arr1;
	@arr2=sort @arr2;
	my $len1=$#arr1;
	my $len2=$#arr2;
	return 0 if $len1!=$len2;
	for(my $i=0;$i<=$len1;$i++){
		return 0 if $arr1[$i] ne $arr2[$i];
	}
	return 1;
}
sub unique{ #用于去重的子程序
	my @array=@_;
	my @tmp;
	push @tmp,shift @array;
	while(@array)
	{
		my $array=shift @array;
		my $same_exist=0;
		foreach (@tmp)
		{
			if($array eq $_)
			{
				$same_exist=1;
			}
		}
		if($same_exist==0)
		{
			push @tmp,$array;
		}
	}
	@tmp;
}