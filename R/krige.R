# Return hidden state estimates
ctmm.state <- function(data,model)
{
  return(kalman(data,model,smooth=TRUE))
}


# Krige at one time or over list of times
krige <- function(Krige, t)
{
  # sampled times
  T <- Krige$t
  N <- length(T)
  dt <- Krige$dt

  # Kriged time indices (bisection lookup)
  n <- length(t)
  I <- array(1,n)
  for(i in 1:n) 
  {
    I.min <- I[i]
    I.max <- N
    
    # iterate two methods:
    # biased/informed bisection that is exact for evenly sampled data
    # naive bisection that always converges
    while(I.max-I.min>1)
    {
      # VANILLA BISECTION
      dI <- floor((I.max-I.min)/2)
      I.test <- I.min + dI
      if(T[I.test]<=t) { I.min <- I.test }
      
      I.test <- I.max - dI
      if(T[I.test]<=t) { I.min <- I.test }
      
      # BIASED/INFORMED BISECTION  
      # How many steps we need to push forward for evenly sampled data
      I.test <- I.min + floor((t[i]-T[I.min])/dt)
      if(T[I.test]<=t) { I.min <- I.test }
      
      # How many steps we need to push backwards for evenly sampled data
      I.test <- I.max - floor((T[I.max]-t[i])/dt)
      if(T[I.test]>=t) { I.max <- I.test }
    }

    # Update max-min
    I[i] <- I.min
    if(i<n) { I[i+1] <- I[i] }
  }
  
  
  
}


# Krige density