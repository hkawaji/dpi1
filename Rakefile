
#---------------------------
#   general setting
#---------------------------
genome = "../in/hg19.chrom.sizes"
ctss_path = "../in/*.bed.gz"
out_path = "../dpiout/"
#---------------------------


#---------------------------
#   parameters
#---------------------------
### for tag clusters
cluster_dist = 20
ctss_max_counts_threshold = 2

### for decomposition
noise_subtraction_ratio = 0.1
length_to_decompose = 50
count_to_decompose = 50
n_comp_upper_bound = 5

### for smoothing, after decomposition
gaussian_window_size_half = 5
#---------------------------


###
### setting parameters from environmental variables (dpi_*)
### ('identify_tss_peaks.sh' override the parameters above)
###
%w(
  dpi_genome dpi_ctss_path dpi_out_path
  dpi_cluster_dist dpi_ctss_max_counts_threshold
  dpi_noise_subtraction_ratio dpi_length_to_decompose
  dpi_count_to_decompose dpi_n_comp_upper_bound
  dpi_gaussian_window_size_half
).each do |k|
  if ENV.key?(k)
    n = k.sub("dpi_","")
    v = ENV[k]
    if ENV.key?("curr_dir")
      if v.match("^/")
        # do nothing, since absolute path
      elsif k.match("path")
        v = ENV["curr_dir"] + "/" + v 
      else
        # do nothing, since 'v' is not path
      end
    end
    eval( "#{n} = '#{v}'" ) 
  end
end





task :prep do |t|
  sh %! mkdir -p #{out_path}/outCounts !
  sh %! mkdir -p #{out_path}/outPooled !
  sh %! mkdir -p #{out_path}/outTpm    !
  sh %! mkdir -p #{out_path}/log       !
end

task :clean do |t|
  sh %! rm -rf #{out_path}/out_pooled/* !
end

task :clean_all do |t|
  sh %! rm -rf #{out_path}/outPooled/* !
  sh %! rm -rf #{out_path}/outCounts/* !
  sh %! rm -rf #{out_path}/outTpm/*    !
  sh %! rm -rf #{out_path}/log/*       !
end


###
### prepare bigWig files for individual samples
###
buf = Hash.new

FileList[ctss_path].each do |infile|

  outfile   = "#{out_path}/outCounts/" + File.basename(infile, ".ctss.bed.gz") + ".ctss.fwd.bw"
  outprefix = "#{out_path}/outCounts/" + File.basename(infile, ".ctss.bed.gz") + ".ctss"
  file outfile do |t|
    sh %! ./scripts/ctssBed2bigWig.sh -i #{infile} -o #{outprefix} -g #{genome}!
  end

  outfileTpm   = "#{out_path}/outTpm/" + File.basename(infile, ".ctss.bed.gz") + ".ctss.fwd.bw"
  outprefixTpm = "#{out_path}/outTpm/" + File.basename(infile, ".ctss.bed.gz") + ".ctss"
  file outfileTpm do |t|
    sh %! ./scripts/ctssBed2TpmBigWig.sh -i #{infile} -o #{outprefixTpm} -g #{genome}!
  end

  task :bw => outfile
  task :bw => outfileTpm
end

###
### pool (CTSS)
###
%w(fwd rev).each do |strand|
  ctssTotalCounts = "#{out_path}/outPooled/ctssTotalCounts.#{strand}.bw"
  ctssMaxCounts   = "#{out_path}/outPooled/ctssMaxCounts.#{strand}.bw"
  ctssMaxTpm      = "#{out_path}/outPooled/ctssMaxTpm.#{strand}.bw"

  file ctssTotalCounts do |t|
    sh %! ulimit -n 2048; bigWigMerge #{out_path}/outCounts/*#{strand}.bw /dev/stdout | sort -k1,1 -k2,2n > #{t.name}.tmp !
    sh %! bedGraphToBigWig #{t.name}.tmp #{genome} #{t.name} !
    sh %! rm -f #{t.name}.tmp !
  end

  file ctssMaxCounts do |t|
    sh %! ulimit -n 2048; bigWigMerge -max #{out_path}/outCounts/*#{strand}.bw /dev/stdout | sort -k1,1 -k2,2n > #{t.name}.tmp !
    sh %! bedGraphToBigWig #{t.name}.tmp #{genome} #{t.name} !
    sh %! rm -f #{t.name}.tmp !
  end

  file ctssMaxTpm do |t|
    sh %! ulimit -n 2048; bigWigMerge -max #{out_path}/outTpm/*#{strand}.bw /dev/stdout | sort -k1,1 -k2,2n > #{t.name}.tmp !
    sh %! bedGraphToBigWig #{t.name}.tmp #{genome} #{t.name} !
    sh %! rm -f #{t.name}.tmp !
  end

  task :pool_ctss => ctssTotalCounts
  task :pool_ctss => ctssMaxCounts
  task :pool_ctss => ctssMaxTpm
end



###
### prepare starting tag cluster
###
file "#{out_path}/outPooled/tc.bed.gz" => :pool_ctss do |t|
  fwdTotal="#{out_path}/outPooled/ctssTotalCounts.fwd.bw"
  revTotal="#{out_path}/outPooled/ctssTotalCounts.rev.bw"

  fwdMax="#{out_path}/outPooled/ctssMaxCounts.fwd.bw"
  revMax="#{out_path}/outPooled/ctssMaxCounts.rev.bw"

  sh %! bigWigToBedGraph #{fwdMax} /dev/stdout \
      | awk --assign ctss_max_counts_threshold=#{ctss_max_counts_threshold} 'BEGIN{OFS="\t"}{if($4 >=ctss_max_counts_threshold){print $1,$2,$3,$4}}' \
      | ./scripts/largeMergeBed.sh -g  #{genome} -d #{cluster_dist} \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4",+",1000,"+"}' \
      | grep -v '_' \
      | bigWigAverageOverBed #{fwdTotal} /dev/stdin  /dev/stdout \
      | awk '{printf "%s\t%i\\n",$1,$4}' \
      > #{t.name}.tmp !

  sh %! bigWigToBedGraph #{revMax} /dev/stdout \
      | awk --assign ctss_max_counts_threshold=#{ctss_max_counts_threshold} 'BEGIN{OFS="\t"}{if($4 >=ctss_max_counts_threshold){print $1,$2,$3,$4}}' \
      | ./scripts/largeMergeBed.sh -g  #{genome} -d #{cluster_dist} \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4",-",1000,"-"}' \
      | grep -v '_' \
      | bigWigAverageOverBed #{revTotal} /dev/stdin  /dev/stdout \
      | awk '{printf "%s\t%i\\n",$1,$4}' \
      >> #{t.name}.tmp !

  sh %! cat #{t.name}.tmp \
      | sed -e 's/[:|,|..]/\t/g' \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1":"$2".."$3","$4,$5,$4}' \
      | sort -k1,1 -k2,2n \
      | gzip -c > #{t.name} !

  sh %! rm -f #{t.name}.tmp !
end
task :tc => "#{out_path}/outPooled/tc.bed.gz"


###
### split tag cluster into the decomposition target (termed "long") and the others (termed "short")
###

### for long
file "#{out_path}/outPooled/tc.long" do |t|
  sh %! rm -rf #{t.name} !
  sh %! mkdir -p #{t.name} !
  sh %! gunzip -c #{out_path}/outPooled/tc.bed.gz \
      | awk '{if( (($3 - $2) > #{length_to_decompose}) && ($5 > #{count_to_decompose}) ){print}}' \
      | gzip -c > #{t.name}.bed.gz !
  sh %! gunzip -c #{t.name}.bed.gz | split --suffix-length=5 --lines=1000  - #{t.name}/ !
end
task :tc_long => "#{out_path}/outPooled/tc.long"

### for short
file "#{out_path}/outPooled/tc.short.bed.gz" do |t|
  sh %! gunzip -c #{out_path}/outPooled/tc.bed.gz \
      | awk '{if( (($3 - $2) <= #{length_to_decompose}) || ($5 <= #{count_to_decompose}) ){print}}' \
      | grep -- '+$' \
      | ./scripts/bed2peakBed9_with_boundaryTrimming.sh -r #{noise_subtraction_ratio} -b #{out_path}/outPooled/ctssTotalCounts.fwd.bw -c '255,0,0' \
      > #{t.name}.tmp !
  sh %! gunzip -c #{out_path}/outPooled/tc.bed.gz \
      | awk '{if( (($3 - $2) <= #{length_to_decompose}) || ($5 <= #{count_to_decompose}) ){print}}' \
      | grep -- '-$' \
      | ./scripts/bed2peakBed9_with_boundaryTrimming.sh -r #{noise_subtraction_ratio} -b #{out_path}/outPooled/ctssTotalCounts.rev.bw -c '0,0,255' \
      >> #{t.name}.tmp !
  sh %! sort -k1,1 -k2,2n #{t.name}.tmp | gzip -c > #{t.name} !
  sh %! rm -f #{t.name}.tmp !
end
task :tc_short => "#{out_path}/outPooled/tc.short.bed.gz"



###
### simple smoothing, followed by threshoulding (without decomposition)
###
file "#{out_path}/outPooled/tc.long.spi.bed.gz" do |t|
  infile_base = "#{out_path}/outPooled/tc.long.bed.gz"
  infile_ctss_pref = "#{out_path}/outPooled/ctssTotalCounts"
  outfile = "#{out_path}/outPooled/tc.long.spi.bed.gz"
  path = "#{out_path}/outPooled"
  pattern="ctssTotalCounts.*.bw$"
  sh %! echo 'source("./scripts/decompose_peakclustering.R");
              main_spi(
                infile_base="#{infile_base}", 
                infile_ctss_pref="#{infile_ctss_pref}",
                gaussian_window_size_half=#{gaussian_window_size_half},
                length_to_decompose = #{length_to_decompose} ,
                noise_subtraction_ratio = #{noise_subtraction_ratio} )' \
      | R --slave \
      | gzip -c > #{t.name} !
end
task :spi_long => "#{out_path}/outPooled/tc.long.spi.bed.gz" 


file "#{out_path}/outPooled/tc.spi_merged.bed.gz" => [:spi_long] do |t|

  sh %! gunzip -c \
          #{out_path}/outPooled/tc.long.spi.bed.gz \
      | grep -v ^# \
      | awk '{if($6=="+"){print}}' \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3, $1":"$2".."$3","$4, 1000,$6 }' \
      | ./scripts/bed2peakBed9.sh -b #{out_path}/outPooled/ctssTotalCounts.fwd.bw -c "255,0,0" \
      > #{t.name}.tmp !

  sh %! gunzip -c \
          #{out_path}/outPooled/tc.long.spi.bed.gz \
      | grep -v ^# \
      | awk '{if($6=="-"){print}}' \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3, $1":"$2".."$3","$4, 1000,$6 }' \
      | ./scripts/bed2peakBed9.sh -b #{out_path}/outPooled/ctssTotalCounts.rev.bw -c "0,0,255" \
      >> #{t.name}.tmp !

  sh %! gunzip -c #{out_path}/outPooled/tc.short.bed.gz >> #{t.name}.tmp !
  sh %! sort -k1,1 -k2,2n #{t.name}.tmp \
      | awk 'BEGIN{OFS="\t"}{print $0,1,$3-$2",",0","}' \
      | gzip -c > #{t.name} !
  sh %! rm -f #{t.name}.tmp !
end
task :spi_merge => "#{out_path}/outPooled/tc.spi_merged.bed.gz"



###
### decomposition, smoothing, and threshoulding
###
subcluster = Hash.new
subcluster[:decompose_smoothing] = []
FileList["#{out_path}/outPooled/tc.long/a*"].each do |infile|
  next if infile.match(/\.gz$/)
  next if infile.match(/\.bed$/)
  infile_base = infile

  outfile = "#{infile}.decompose_smoothing.bed"
  path = "#{out_path}/outCounts"
  pattern=".bw$"
  file outfile do |t|
    fh = open('| mktemp -p "${TMPDIR:-/tmp}" -d dpi.XXXXXXXX ')
    tmpdir = fh.gets.chomp
    fh.close()
    local_path="#{tmpdir}/#{File.basename(path)}"
    sh %! trap '[[ #{tmpdir} ]] && rm -rf #{tmpdir}' 0 1 2 3 15 ; \
          cp -rp #{path} #{tmpdir}/ ; \
          echo 'source("./scripts/decompose_peakclustering.R");
                main_dpi(
                  infile_base="#{infile_base}", path="#{local_path}", pattern="#{pattern}",
                  exclude_prefix = "xxxxxxxxxxxx\.",
                  gaussian_window_size_half=#{gaussian_window_size_half},
                  n.comp.upper_bound = #{n_comp_upper_bound},
                  length_to_decompose = #{length_to_decompose} ,
                  noise_subtraction_ratio = #{noise_subtraction_ratio} )' \
        | R --slave \
        > #{t.name} !
  end
  task :decompose_smoothing => outfile
  subcluster[:decompose_smoothing] << outfile
end

subcluster.keys.each do |key|

  outfile = "#{out_path}/outPooled/tc.long.#{key.to_s}.bed.gz"

  file outfile do |t|
    colors = {
      "+" => "255,0,0",
      "-" => "0,0,255",
    }
    command = %! cat #{subcluster[key].join(' ')} \
               | grep -v '^#' \
               | awk '{if($1 \!~ ":"){print}}' \
               | awk '{printf "%s\t%i\t%i\t%s\t%i\t%s\\n",$1,$2,$3,$4,$5,$6}' \
               | sort -k 1,1 -k 2,2n !
    open("| gzip -c > #{t.name}", "w") do |ofh|
    open("| #{command}")  do |ifh|
      ifh.each do |line|
        #chr1    566854  566953  chr1:566727..567119,+;peak:chr1:566874..566875,+;peakCounts:348 1000    +
        #chr1    909815  910037  chr1:909815..910037,+   1000    +
        cols = line.chomp.split(/\t/)
        if m = cols[3].match(/peak:chr.*:(\d+)\.\.(\d+),/)
          m = m.to_a
          cols << m[1]
          cols << m[2]
        else
          cols << cols[1]
          cols << cols[2]
        end
        cols << colors[cols[5]]
        ofh.puts cols.join("\t")
        #puts cols.join("\t")
      end # line
    end # ifh
    end # ofh
  end

  (1 .. n_comp_upper_bound).each do |n|
    %w( fwd rev ).each do |strand|
      strand_symbol = { "fwd" => "+", "rev" => "-"}
      outfile3 = "#{out_path}/outPooled/tc.long.#{key.to_s}.component#{n}_ctss.#{strand}.bedGraph.gz"
      file outfile3 do |t|
        sh %! cat #{subcluster[key].join(' ')} \
            | grep -v '^#' \
            | awk --assign n=#{n} '{if(($1 ~ ":") && (NF > n)){print}}' \
            | cut -f 1,#{n + 1} \
            | sed 's/:/\t/;s/,/\t/;s/\\.\\./\t/' \
            | awk --assign strand=#{strand_symbol[strand]} '{if (($4 == strand) && ($5 \!= "NA")){printf "%s\t%i\t%i\t%i\\n", $1,$2,$3,$5} }' \
            | sort -k 1,1 -k 2,2n  \
            | gzip -c > #{t.name} !
      end
      task :comp_ctss => outfile3
    end
  end

end



subcluster.keys.each do |key|
  infile = "#{out_path}/outPooled/tc.long.#{key.to_s}.bed.gz"
  outfile = "#{out_path}/outPooled/tc.long.#{key.to_s}.merged.bed.gz"
  file outfile => infile do |t|

    sh %! gunzip -c #{infile} \
        | awk '{if($6=="+"){print}}' \
        | cut -f 1-6 \
        | mergeBed -i stdin \
        | awk 'BEGIN{OFS="\t";strand="+"}{print $1,$2,$3, $1":"$2".."$3","strand, 1000,strand }' \
        | ./scripts/bed2peakBed9.sh -b #{out_path}/outPooled/ctssTotalCounts.fwd.bw -c "255,0,0" \
        > #{t.name}.tmp !

    sh %! gunzip -c #{infile} \
        | awk '{if($6=="-"){print}}' \
        | cut -f 1-6 \
        | mergeBed -i stdin \
        | awk 'BEGIN{OFS="\t";strand="-"}{print $1,$2,$3, $1":"$2".."$3","strand, 1000,strand }' \
        | ./scripts/bed2peakBed9.sh -b #{out_path}/outPooled/ctssTotalCounts.rev.bw -c "0,0,255" \
        >> #{t.name}.tmp !

    sh %! sort -k1,1 -k2,2n #{t.name}.tmp | gzip -c > #{t.name} !
    sh %! rm -f #{t.name}.tmp !

  end

  outfile_final = "#{out_path}/outPooled/tc.#{key.to_s}_merged.bed.gz"
  file outfile_final => [outfile] do |t|
    sh %! gunzip -c #{outfile} #{out_path}/outPooled/tc.short.bed.gz \
        | sort -k1,1 -k2,2n -k6,6 \
        | awk 'BEGIN{OFS="\t"}{print $0,1,$3-$2",",0","}' \
        | gzip -c > #{t.name} !
  end
  task :final => outfile_final
end


### permissive
[3].each do |cutoff|
  outfile_dpi = "#{out_path}/outPooled/tc.decompose_smoothing_merged.ctssMaxCounts#{cutoff}.bed.gz"
  file outfile_dpi do |t|
    sh %! bigWigToBedGraph #{out_path}/outPooled/ctssMaxCounts.fwd.bw /dev/stdout \
        | awk --assign cutoff=#{cutoff} 'BEGIN{OFS="\t"}{if($4 >= cutoff ){print $1,$2,$3,".",$4,"+"}}' > #{t.name}.tmp !
    sh %! bigWigToBedGraph #{out_path}/outPooled/ctssMaxCounts.rev.bw /dev/stdout \
        | awk --assign cutoff=#{cutoff} 'BEGIN{OFS="\t"}{if($4 >= cutoff ){print $1,$2,$3,".",$4,"-"}}' >> #{t.name}.tmp !
    sh %! gunzip -c #{out_path}/outPooled/tc.decompose_smoothing_merged.bed.gz \
        | intersectBed -s -wa -u -a stdin -b #{t.name}.tmp \
        | gzip -c > #{t.name} !
    sh %! rm -f #{t.name}.tmp !
  end
  task :threshold => outfile_dpi


  outfile_spi = "#{out_path}/outPooled/tc.spi_merged.ctssMaxCounts#{cutoff}.bed.gz"
  file outfile_spi => :spi_merge do |t|
    sh %! bigWigToBedGraph #{out_path}/outPooled/ctssMaxCounts.fwd.bw /dev/stdout \
        | awk --assign cutoff=#{cutoff} 'BEGIN{OFS="\t"}{if($4 >= cutoff ){print $1,$2,$3,".",$4,"+"}}' > #{t.name}.tmp !
    sh %! bigWigToBedGraph #{out_path}/outPooled/ctssMaxCounts.rev.bw /dev/stdout \
        | awk --assign cutoff=#{cutoff} 'BEGIN{OFS="\t"}{if($4 >= cutoff ){print $1,$2,$3,".",$4,"-"}}' >> #{t.name}.tmp !
    sh %! gunzip -c #{out_path}/outPooled/tc.spi_merged.bed.gz \
        | intersectBed -s -wa -u -a stdin -b #{t.name}.tmp \
        | gzip -c > #{t.name} !
    sh %! rm -f #{t.name}.tmp !
  end
  task :spi => outfile_spi

end

### robust
[[11, 1]].each do |cutoff,cutoff_tpm|
  outfile_dpi = "#{out_path}/outPooled/tc.decompose_smoothing_merged.ctssMaxCounts#{cutoff}_ctssMaxTpm#{cutoff_tpm}.bed.gz"
  file outfile_dpi do |t|
    sh %! bigWigToBedGraph #{out_path}/outPooled/ctssMaxCounts.fwd.bw /dev/stdout \
        | awk --assign cutoff=#{cutoff} 'BEGIN{OFS="\t"}{if($4 >= cutoff ){print $1,$2,$3,".",$4,"+"}}' > #{t.name}.tmp !
    sh %! bigWigToBedGraph #{out_path}/outPooled/ctssMaxCounts.rev.bw /dev/stdout \
        | awk --assign cutoff=#{cutoff} 'BEGIN{OFS="\t"}{if($4 >= cutoff ){print $1,$2,$3,".",$4,"-"}}' >> #{t.name}.tmp !
    sh %! gunzip -c #{out_path}/outPooled/tc.decompose_smoothing_merged.bed.gz \
        | intersectBed -s -wa -u -a stdin -b #{t.name}.tmp \
        | gzip -c > #{t.name}.tmp2.gz !

    sh %! bigWigToBedGraph #{out_path}/outPooled/ctssMaxTpm.fwd.bw /dev/stdout \
        | awk --assign cutoff_tpm=#{cutoff_tpm} 'BEGIN{OFS="\t"}{if($4 >= cutoff_tpm ){print $1,$2,$3,".",$4,"+"}}' > #{t.name}.tmp !
    sh %! bigWigToBedGraph #{out_path}/outPooled/ctssMaxTpm.rev.bw /dev/stdout \
        | awk --assign cutoff_tpm=#{cutoff_tpm} 'BEGIN{OFS="\t"}{if($4 >= cutoff_tpm ){print $1,$2,$3,".",$4,"-"}}' >> #{t.name}.tmp !
    sh %! intersectBed -s -wa -u -a #{t.name}.tmp2.gz  -b #{t.name}.tmp \
        | gzip -c > #{t.name} !

    sh %! rm -f #{t.name}.tmp #{t.name}.tmp2.gz !
  end
  task :threshold => outfile_dpi


  outfile_spi = "#{out_path}/outPooled/tc.spi_merged.ctssMaxCounts#{cutoff}_ctssMaxTpm#{cutoff_tpm}.bed.gz"
  file outfile_spi => :spi_merge do |t|
    sh %! bigWigToBedGraph #{out_path}/outPooled/ctssMaxCounts.fwd.bw /dev/stdout \
        | awk --assign cutoff=#{cutoff} 'BEGIN{OFS="\t"}{if($4 >= cutoff ){print $1,$2,$3,".",$4,"+"}}' > #{t.name}.tmp !
    sh %! bigWigToBedGraph #{out_path}/outPooled/ctssMaxCounts.rev.bw /dev/stdout \
        | awk --assign cutoff=#{cutoff} 'BEGIN{OFS="\t"}{if($4 >= cutoff ){print $1,$2,$3,".",$4,"-"}}' >> #{t.name}.tmp !
    sh %! gunzip -c #{out_path}/outPooled/tc.spi_merged.bed.gz \
        | intersectBed -s -wa -u -a stdin -b #{t.name}.tmp \
        | gzip -c > #{t.name}.tmp2.gz !

    sh %! bigWigToBedGraph #{out_path}/outPooled/ctssMaxTpm.fwd.bw /dev/stdout \
        | awk --assign cutoff_tpm=#{cutoff_tpm} 'BEGIN{OFS="\t"}{if($4 >= cutoff_tpm ){print $1,$2,$3,".",$4,"+"}}' > #{t.name}.tmp !
    sh %! bigWigToBedGraph #{out_path}/outPooled/ctssMaxTpm.rev.bw /dev/stdout \
        | awk --assign cutoff_tpm=#{cutoff_tpm} 'BEGIN{OFS="\t"}{if($4 >= cutoff_tpm ){print $1,$2,$3,".",$4,"-"}}' >> #{t.name}.tmp !
    sh %! intersectBed -s -wa -u -a #{t.name}.tmp2.gz  -b #{t.name}.tmp \
        | gzip -c > #{t.name} !

    sh %! rm -f #{t.name}.tmp #{t.name}.tmp2.gz !
  end
  task :spi => outfile_spi


end



