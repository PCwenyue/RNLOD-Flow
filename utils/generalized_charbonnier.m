function y = generalized_charbonnier(x, sigma, type)
%CHARBONNIER  GENERALIZED_Charbonnier robust function.
%   CHARBONNIER(X, SIGMA, TYPE) evaluates the Charbonnier robust function
%   with sigma SIGMA at point(s) X.  
%   TYPE selects the evaluation type:
%    - 0: function value
%    - 1: first derivative
%    - 2: second derivative
%  
%   This is a private member function of the class 'robust_function'. 
%
% Authors: Deqing Sun, Department of Computer Science, Brown University
%          Stefan Roth, Department of Computer Science, TU Darmstadt
% Contact: dqsun@cs.brown.edu, sroth@cs.tu-darmstadt.de
% $Date: $
% $Revision: $

% Copyright 2004-2007, Brown University, Providence, RI. USA
% Copyright 2007-2010 TU Darmstadt, Darmstadt, Germany.
% 
%                          All Rights Reserved
% 
% All commercial use of this software, whether direct or indirect, is
% strictly prohibited including, without limitation, incorporation into in
% a commercial product, use in a commercial service, or production of other
% artifacts for commercial purposes.     
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for research purposes is hereby granted without fee,
% provided that the above copyright notice appears in all copies and that
% both that copyright notice and this permission notice appear in
% supporting documentation, and that the name of the author and Brown
% University not be used in advertising or publicity pertaining to
% distribution of the software without specific, written prior permission.        
%
% For commercial uses contact the Technology Venture Office of Brown University
% 
% THE AUTHOR AND BROWN UNIVERSITY DISCLAIM ALL WARRANTIES WITH REGARD TO
% THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
% FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR OR
% BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
% DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
% PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
% ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
% THIS SOFTWARE.        

  sig  = sigma(1);
  a    = sigma(2);

  switch (type)
   case 0
    y = (sig^2 + x.^2).^a;
   case 1
    y = 2*a*x.*(sig^2 + x.^2).^(a-1);
    %y = x ./ sqrt(1 + x.^2 / sigma^2);
   case 2
    y = 2*a*(sig^2 + x.^2).^(a-1);
     case 3
    y = 2*a*(sig^2 + x).^(a-1);
  end