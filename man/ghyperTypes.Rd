\name{ghyper.types}
\alias{ghyper.types}


\title{Kemp and Kemp generalized hypergeometric types}
\description{
	Generalized hypergeometric types as given by Kemp and Kemp	
}

\section{Two-way table}{

The basic representation is in terms of a two-way table:

\tabular{ccc}{
x \tab  k-x \tab  k \cr
a-x \tab  b-k+x \tab  N-k \cr
  a \tab     b \tab  N
}

and the associated hypergeometric probability \eqn{P(x)=C_x^a C_{k-x}^b / C_k^N}{P(x)=choose(a, x) choose(b, k-x) / choose(N, k)}.

The types are classified according to ranges of a, k, and N.

}
\section{Kemp and Kemp types}{

Minor modifications in the definition of three of the types have been made to avoid numerical difficulties. Note, J denotes a nonnegative integer.


\tabular{ll}{
[Classic] \tab \cr
	\tab \eqn{0<a, 0<N, 0<k} \cr
 	\tab	integers: a, N, k.\cr
	\tab	\eqn{max(0,a+k-N)  \le  x  \le  min(a,k)}{max(0,a+k-N) <= x <= min(a,k)}\cr
[IA(i)] (Real classic) \tab at least one noninteger parameter\cr
	\tab \eqn{0<a, 0<N, 0<k, k-1<a<N-(k-1)} \cr
	\tab	integer: k \cr
	\tab \eqn{	0  \le  x  \le  a}{0 <= x <= k}\cr

[IA(ii)] (Real classic)\tab at least one noninteger parameter\cr
	\tab \eqn{0<a, 0<N, 0<k, a-1<k<N-(a-1)} \cr
	\tab integer: a \cr 
	\tab \eqn{0  \le  x  \le  a}{0 <= x <= a} \cr 
	\tab Interchanging a and k transforms this to type IA(i) \cr

[IB] \tab \cr
	\tab \eqn{0<a, 0<N, 0<k, a+k-1<N, J < (a,k) < J+1} \cr
	\tab integer: \eqn{0 \le J}{0 <= J} \cr
	\tab non-integer: a, k \cr
	\tab \eqn{0  <=  x  \dots }{0 <= x \ldots}\cr
	
	\tab NOTE: Kemp and Kemp specify \eqn{-1<N}{-1<N}. \cr
	\tab No practical applications for this distribution. \cr
[IIA] (negative hypergeometric) \tab \cr
	\tab \eqn{a<0, N<a-1,0<k} \cr
	\tab integer: k \cr
	\tab \eqn{0  \le  x  \le  k}{0 <= x <= k} \cr
	
	\tab NOTE: Kemp and Kemp specify \eqn{N<a,  N \ne a-1}{N<a and  N!=(a-1)} \cr
[IIB] \tab \cr
	\tab \eqn{a<0, -1<N<k+a-1, 0<k, J < (k,k+a-1-N) < J+1} \cr

	\tab non-integer: k \cr
	\tab integer: \eqn{0 \le J}{0 <= J} \cr
	\tab \eqn{0  \le  x  ....}{0 <= x \ldots} \cr

	\tab This is a very strange distribution.  Special calculations were used.\cr

	\tab Note: No practical applications.\cr
[IIIA] (negative hypergeometric) \tab \cr
	\tab \eqn{0<a,N<k-1,k<0}\cr
	\tab integer: a\cr
	\tab \eqn{0 \le x \le a}{0  <=  x  <=  a}\cr
	
	\tab Interchanging a and k transforms this to type IIA \cr

	\tab NOTE: Kemp and Kemp specify \eqn{N<k, N \ne k-1}{N<k and N != k-1} \cr

[IIIB]\tab \cr
	\tab  \eqn{0<a,-1<N<a+k-1,k<0, J<(a,a+k-1-N)<J+1} \cr
	\tab non integer: a\cr
	\tab integer: \eqn{0 \le J}{0 <= J}\cr
	\tab \eqn{0 \le x \dots }{0  <=  x \ldots}\cr

	\tab Interchanging a and k transforms this to type IIB \cr
	\tab Note: No practical applications\cr

[IV] (Generalized Waring) \tab \cr
	\tab \eqn{a<0,-1<N, k<0}\cr
	\tab \eqn{0 \le x \dots}{0  <=  x  \ldots}

}


}


\references{

 
Kemp, C.D., and Kemp, A.W. (1956). Generalized hypergeometric distributions. \emph{Jour. Roy. Statist. Soc. B.} \bold{18.} 202-211.
\bold{39.}  887-895.
}

\author{
	Bob Wheeler \email{bwheelerg@gmail.com}
}
\keyword{distribution}