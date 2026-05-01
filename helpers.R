get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  
  # Check if run via Rscript
  match <- grep(file_arg, args)
  if (length(match) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", args[match]))))
  }
  
  # Check if run in RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }
  
  # Fallback for source()
  return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
}

get_lst_depth = function(this, this_depth=0) {
  if (is.data.frame(this) | !(is.list(this))) {
    return(this_depth)
  } else {
    return(max(unlist(lapply(this, get_lst_depth, this_depth=this_depth+1))))    
  }
}

apply_last_level = function(lst, f) {
  # Recursion
  apply_last_level_helper = function(this_lst, this_names) {
    if (is.list(this_lst) & !(is.data.frame(this_lst))) {
      res = lapply(seq_along(this_lst), function(i) {
        return(apply_last_level_helper(this_lst[[i]], names(this_lst)[i]))
      })
      names(res) = names(this_lst)
      return(res)
    } else {
      return(f(this_lst, this_names))
    }
  }
  res_lst = apply_last_level_helper(lst, c())
  names(res_lst) = names(lst)
  return(res_lst)
}

apply_last_level_together = function(lst1, lst2, f) {
  apply_last_level_together_helper = function(this_lst1, this_lst2) {
    if (is.list(this_lst1) & 
        !(is.data.frame(this_lst1)))  {
      if (is.list(this_lst2) & 
          !(is.data.frame(this_lst2))) {
        if (all(names(this_lst1) == names(this_lst2))) {
          # Step down both nested lists
          res = mapply(function(a, b) {
            return(apply_last_level_together_helper(a, b))
          }, a=this_lst1, b=this_lst2, SIMPLIFY=FALSE)
          names(res) = names(this_lst1)
          return(res)
        } else {
          # Step down first nested list only
          res = lapply(this_lst1, function(x) {
            return(apply_last_level_together_helper(x, this_lst2))
          })
          names(res) = names(this_lst1)
          return(res)
        }
      } else {
        # Step down first nested list only
        res = lapply(this_lst1, function(x) {
          return(apply_last_level_together_helper(x, this_lst2))
        })
        names(res) = names(this_lst1)
        return(res)
      }
    } else {
      return(f(this_lst1, this_lst2))
    }
  }
  return(apply_last_level_together_helper(lst1, lst2))
}

read_table = function(tbl_path, sep="\t", header=TRUE, skip=0) {
  return(
    read.table(tbl_path, sep=sep, header=header, 
               stringsAsFactors=FALSE, quote="", comment.char="",
               skip=skip)
  )
}

write_table = function(tbl, tbl_path, sep="\t", header=TRUE) {
  write.table(tbl, tbl_path,
              sep=sep, row.names=FALSE, quote=FALSE, col.names=header)
}

rm_null_in_list = function(lst) {
  return(lst[!(sapply(lst, is.null))])
}

rm_ensembl_version = function(ensembl_ids) {
  return(gsub("^([^\\.]*)\\..*$", "\\1", ensembl_ids))
}

rm_chr = function(chr_names) {
  return(gsub("^chr(.*)", "\\1", chr_names))
}
