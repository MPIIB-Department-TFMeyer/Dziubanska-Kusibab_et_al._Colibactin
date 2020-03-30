
all_pwm_raw = list()

dreme_files = list.files(path=dreme_result_folder, pattern="dreme.txt", recursive = T, full.names = T)

for (ff in dreme_files) {
  tmp = unlist(strsplit(ff,"/"))
  s = tmp[length(tmp)-1]
  
  if (! s %in% dreme_libs) next
  
  tt = scan(file=ff, sep="\n", blank.lines.skip = F, what=character(), quiet = T)
  startlines = grep("^letter-probability matrix:", tt)
  blank_lines = grep("^$", tt)
  endlines = unlist(Map(function(x) {blank_lines[min(which(blank_lines > x))]}, startlines))
  
  blocks = mapply(function(s, e) paste(tt[s:e], collapse="\n"), startlines, endlines-1, SIMPLIFY = F)
  names(blocks) = paste0("m", 1:length(blocks))
  all_pwm_raw[[s]] = blocks
}


require(Biostrings)

all_pwm = list()
for (s in names(all_pwm_raw)) {
  all_pwm[[s]] = list()
  for (mn in names(all_pwm_raw[[s]])) {
    raw = all_pwm_raw[[s]][[mn]]
    mat = read.table(text = raw, sep=" ", skip=1)
    colnames(mat) = c("A","C","G","T")
    par_raw = scan(text=gsub("letter-probability matrix: ", "", raw), what=character(), nlines=1, sep=" ", quiet=T)
    par_vec = par_raw[seq(2, length(par_raw), 2)]
    names(par_vec) = gsub("=", "", par_raw[seq(1, length(par_raw), 2)])
    total_sites = as.integer(par_vec["nsites"])
    mat_tmp = round(t(mat) * total_sites)
    mode(mat_tmp) = "integer"
    col_sums = apply(mat_tmp, 2, sum)
    col_missing = total_sites - col_sums
    if (!all(col_missing == 0)) {
      # simple and dirty: add the missing to the row with the highest value
      mat_tmp = mapply(function(column_ind, diff) {v = mat_tmp[, column_ind]; if(diff >= 0){ri = which(max(v)==v); v[ri] = v[ri]+diff}; return(v)}, 1:ncol(mat_tmp), col_missing)
      #stop("column sum differences detected");
    } 
    all_pwm[[s]][[mn]] = list("pwm"=PWM(mat_tmp), "rel"=t(mat), "count"=mat_tmp, "parameters"=par_vec)
  }
}
