library(data.table)
library(ggplot2)
max.angle <- 360
check_angle <- function(x){
  stopifnot(0 <= x & x < max.angle)
}
geodesic_dist <- function(data.vec, param.vec){
  check_angle(data.vec)
  check_angle(param.vec)
  d <- abs(data.vec-param.vec)
  ifelse(d>max.angle/2, max.angle-d, d)
}
bkp <- function(Linear_diff, Constant_diff, param){
  data.table(param, Linear_diff, Constant_diff)
}
bkps_for_data <- function(data.val){
  if(data.val == 0){
    bkp(-2, max.angle, max.angle/2)
  }else if(data.val < max.angle/2)rbind(
    bkp(2, -2*data.val, data.val),
    bkp(-2, max.angle+2*data.val, data.val+max.angle/2)
  )else if(data.val == max.angle/2){
    bkp(2, -max.angle, max.angle/2)
  }else rbind(
    bkp(-2, 2*data.val-max.angle, data.val-max.angle/2),
    bkp(2, -2*data.val, data.val)
  )
}
piece <- function(Linear, Constant, min_param, max_param){
  data.table(Linear, Constant, min_param, max_param)
}
pieces_for_data <- function(data.val){
  if(data.val == 0)rbind(
    piece(1, 0, 0, max.angle/2),
    piece(-1, max.angle, max.angle/2, max.angle)
  )else if(data.val < max.angle/2)rbind(
    piece(-1, data.val, 0, data.val),
    piece(1, -data.val, data.val, data.val+max.angle/2),
    piece(-1, (max.angle+data.val), data.val+max.angle/2, max.angle)
  )else if(data.val == max.angle/2)rbind(
    piece(-1, max.angle/2, 0, max.angle/2),
    piece(1, -max.angle/2, max.angle/2, max.angle)
  )else rbind(
    piece(1, max.angle-data.val, 0, data.val-max.angle/2),
    piece(-1, data.val, data.val-max.angle/2, data.val),
    piece(1, -data.val, data.val, max.angle)
  )
}
pieces_for_data(5)

grid.dt <- data.table(degrees=seq(0,359))
data.list <- list(
  ##start_zero=c(0,50,5,55,200),#TODO
  move_before=c(30,103,230,250,200),
  simple=c(270,180,200,0,90))
loss.dt.list <- list()
diff.dt.list <- list()
ptr.dt.list <- list()
for(Data in names(data.list)){
  data.vec <- data.list[[Data]]
  current.diffs <- data.table()
  pointer_i <- 1
  pointer_param <- max.angle
  pointer_Linear <- 0
  pointer_Constant <- 0
  for(data.i in seq_along(data.vec)){
    data.value <- data.vec[[data.i]]
    right_i <- if(pointer_i>nrow(current.diffs))1 else pointer_i+1
    diff_row <- if(right_i==1){
      data.table(Linear_diff=0, Constant_diff=-max.angle)
    }else current.diffs[pointer_i]
    right_Linear <- pointer_Linear+diff_row$Linear_diff
    is.flat <- right_Linear == 0
    if(nrow(current.diffs) && is.flat){
      print(diff_row)
      right_param <- current.diffs[right_i,param]+if(right_i==1)max.angle else 0
      maybe_over <- (pointer_param+right_param)/2
      mid_param <- if(maybe_over>=max.angle)maybe_over-max.angle else maybe_over
      max_param_over <- data.value+max.angle/2
      max_param <- if(max_param_over>=max.angle)max_param_over-max.angle else max_param_over
      move.right <- if(pointer_param<mid_param){
        pointer_param < max_param & max_param < mid_param
      }else{
        stop("TODO")
      }
      if(move.right){
        print(data.table(pointer_i,pointer_param,pointer_Linear,pointer_Constant,loss=pointer_param*pointer_Linear+pointer_Constant))
        pointer_i <- right_i
        pointer_param <- current.diffs[pointer_i, param]
        pointer_Linear <- right_Linear
        pointer_Constant <- pointer_Constant+diff_row$Constant_diff
        print(data.table(pointer_i,pointer_param,pointer_Linear,pointer_Constant,loss=pointer_param*pointer_Linear+pointer_Constant))
      }
    }
    current.diffs <- rbind(
      current.diffs,
      bkps_for_data(data.value)
    )[order(param)]
    current.diffs <- dcast(
      current.diffs,
      param ~ ., sum,
      value.var=c("Linear_diff","Constant_diff")
    )[!(Linear_diff==0 & Constant_diff==0)]
    diff.dt.list[[paste(Data,data.i)]] <- data.table(
      Data,data.i, current.diffs)
    pointer_i <- if(pointer_param==max.angle){
      nrow(current.diffs)+1
    }else{
      which(current.diffs$param==pointer_param)
    }
    piece_to_add <- pieces_for_data(data.value)[
      min_param <= pointer_param & pointer_param <= max_param]
    pointer_Linear <- pointer_Linear + piece_to_add$Linear
    pointer_Constant <- pointer_Constant + piece_to_add$Constant
    moves <- 0
    ptr.dt.list[[paste(Data,data.i, "before_min")]] <- data.table(
      Data,data.i, Pointer="before", param=pointer_param, moves,
      loss=pointer_param*pointer_Linear+pointer_Constant)
    pointer_at_min <- FALSE
    while(!pointer_at_min){
      right_i <- if(pointer_i>nrow(current.diffs))1 else pointer_i+1
      left_i <- if(pointer_i==1)nrow(current.diffs)+1 else pointer_i-1
      if(pointer_Linear>=0){
        ## increasing, want to move left.
        diff_row <- if(left_i>nrow(current.diffs)){
          data.table(Linear_diff=0, Constant_diff=-max.angle)
        }else current.diffs[left_i]
        moves <- moves+1
        pointer_i <- left_i
        pointer_param <- current.diffs[pointer_i, param]
        pointer_Linear <- pointer_Linear-diff_row$Linear_diff
        pointer_Constant <- pointer_Constant-diff_row$Constant_diff
      }else if(pointer_Linear<0){
        ## decreasing, want to move right.
        diff_row <- if(right_i==1){
          data.table(Linear_diff=0, Constant_diff=-max.angle)
        }else current.diffs[pointer_i]
        right_Linear <- pointer_Linear+diff_row$Linear_diff
        right2_i <- if(right_i>nrow(current.diffs))1 else right_i+1
        diff_right <- if(right2_i==1){
          data.table(Linear_diff=0, Constant_diff=-max.angle)
        }else current.diffs[right_i]
        right2_Linear <- right_Linear+diff_right$Linear_diff
        if(right_Linear<0){#move right.
          moves <- moves+1
          pointer_i <- right_i
          pointer_param <- current.diffs[pointer_i, param]
          pointer_Linear <- right_Linear
          pointer_Constant <- pointer_Constant+diff_row$Constant_diff
        }else if(right_Linear==0 && right2_Linear<0){
          moves <- moves+2
          print(moves)
          pointer_i <- right2_i
          pointer_param <- current.diffs[pointer_i, param]
          pointer_Linear <- right2_Linear
          pointer_Constant <- pointer_Constant+
            diff_row$Constant_diff+diff_right$Constant_diff
        }else{
          pointer_at_min <- TRUE
        }
      }
    }
    after_min <- data.table(
      Data, data.i, Pointer="after", param=pointer_param, moves,
      loss=pointer_param*pointer_Linear+pointer_Constant)
    ptr.dt.list[[paste(Data,data.i, "after_min")]] <- after_min
    ptr.dt.list[[paste(Data,data.i, "after_max")]] <- data.table(
      after_min
    )[, `:=`(
      Pointer="after_max",
      param=ifelse(param<max.angle/2,param+max.angle/2, param-max.angle/2),
      loss=max.angle*data.i/2-loss
    )][]
    loss.dt.list[[paste(Data,data.i)]] <- grid.dt[, data.table(
      Data,data.i,
      data.value,
      param=degrees,
      loss=geodesic_dist(data.value,degrees))]
  }
}

(diff.dt <- rbindlist(diff.dt.list))
(loss.dt <- rbindlist(
  loss.dt.list
)[
, cum.loss := cumsum(loss)
, by=.(Data,param)][])
ptr.dt <- rbindlist(ptr.dt.list)
gg <- ggplot()+
  theme_bw()+
  geom_vline(aes(
    xintercept=param),
    data=diff.dt,
    color="grey")+
  geom_line(aes(
    param, cum.loss),
    data=loss.dt)+
  facet_wrap(~Data+data.i, scales="free", labeller=label_both, ncol=5)+
  geom_point(aes(
    param, loss, fill=Pointer, size=Pointer),
    shape=21,
    data=ptr.dt)+
  scale_fill_manual(values=c(
    before="white",
    after="black",
    after_max="red"))+
  scale_size_manual(values=c(
    before=2,
    after=1,
    after_max=1))+
  scale_x_continuous(
    "parameter",
    breaks=seq(0, 360, by=90))+
  scale_y_continuous(
    "Geodesic loss")
png("figure-pointer-moves.png", width=10, height=4, units="in", res=200)
print(gg)
dev.off()
