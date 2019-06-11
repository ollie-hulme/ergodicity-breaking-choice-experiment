function p = sanePColor(varargin)
%SANEPCOLOR  simple wrapper for pcolor
%
% Unlike the built-in pcolor command, this function does not "cut off" the
% last row and column of the input matrix.  In this way, sanePColor is
% intended to be as easy to use as imagesc, but allows the user to specify
% the x and y coordinates of each cell if desired.  This function is also
% useful as an alternative means of generating images to print to PDF that
% are compatible with OS X's "Preview" PDF viewer (imagesc images appear
% "blurred" when printing to a PDF as a vector graphic and viewed using
% Preview).
%
% NOTE: The imagesc function assumes that each entry in a matrix gets the
% corresponding coordinate in 2-d space.  For example, entry (2,3) in the
% matrix is assigned to the coordinate with x = 3, y = 2 when using
% imagesc.  The pcolor function assumes that entries correspond to *edges*.
% So entry (2,3) in a matrix corresponds to the value between 2 and 3
% (along the x-axis) and between 1 and 2 (along the y-axis).  This is why
% one row and one column are cut off when using pcolor.  sanePColor
% behaves like imagesc (i.e. does not cut off data), but uses the "edge
% assignment" data representation required by pcolor.  sanePColor uses
% linear or logarithmic interpolation to infer the edges automatically.
%
% Usage: p = sanePColor([x,y],z,[logx],[logy]);
%
%INPUTS:
%
%    x: an array of sorted x values.  can also specify a min and max x value.
%       these values correspond to columns of z. [IF THIS ARGUMENT IS USED,
%       MUST ALSO SPECIFY Y VALUES.]
% 
%    y: an array of sorted y values.  can also specify a min and max y value.
%       these values correspond to rows of z.  [IF THIS ARGUMENT IS USED,
%       MUST ALSO SPECIFY X VALUES.]
% 
%    z: a 2d matrix of values.  this matrix determines the color at each
%       point.
% 
% logx: if this optional argument is set to true, the x-axis will plotted
%       in log scale (similar to semilogx).
%
% logy: if this optional argument is set to true, the y-axis will plotted
%       in log scale (similar to semilogy).
%
%OUTPUTS:
%
%    p: a handle to the resulting pcolor image.
%
% EXAMPLE:
%
%   m = membrane;
%   p = sanePColor(m);
%
% SEE ALSO: PCOLOR, IMAGE, IMAGESC, SEMILOGX, SEMILOGY, LOGLOG, PADARRAY
%
%   AUTHOR: JEREMY R. MANNING
%  CONTACT: jeremy@dartmouth.edu
%CHANGELOG
%3-16-10    JRM      Wrote it.
%3-12-12    JRM      Support a more diverse range of input configurations.
%9-21-12    JRM      Use linear and logistic interpolation to estimate data
%                    coordinates more accurately.
%10-6-16    JRM      Support non-linear and non-log axes (Credit: Benjamin
%                    Strom suggestion via Mathworks FileExchange)
%parse arguments
if length(varargin) == 1 %just z
    z = varargin{1};
    x = 1:size(z,2);
    y = 1:size(z,1);
    [logx,logy] = deal(false);    
elseif (length(varargin) >= 4) %x, y, z, logx, and possibly logy
    x = varargin{1};    
    y = varargin{2};    
    z = varargin{3};    
    logx = varargin{4};
    if length(varargin) >= 5, logy = varargin{5}; else logy = false; end
elseif length(varargin) == 2 %z and logx
    z = varargin{1};
    logx = varargin{2};
    assert(islogical(logx),'logx must be a logical');
    if logx
        x = logspace(log10(1),log10(size(z,2)),size(z,2));
    else
        x = 1:size(z,2);
    end
    logy = false;
    y = 1:size(z,1);
else %length(varargin) == 3
    if isempty(varargin)
        fprintf('\nUsage: p = sanePColor([x,y],z,[logx],[logy]);\n');
        fprintf('Type ''help %s'' for more info.\n\n',mfilename);
        p = [];
        return;
    end
    %posibility 1: x, y, z
    if length(varargin{1}) > 1 && length(varargin{2}) > 1
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        [logx,logy] = deal(false);
    %posibility 2: z, logx, and logy
    else
        z = varargin{1};
        logx = varargin{2};
        assert(islogical(logx),'logx must be a logical');
        if logx
            x = logspace(log10(1),log10(size(z,2)),size(z,2));
        else
            x = 1:size(z,2);
        end
        
        logy = varargin{3};
        assert(islogical(logy),'logy must be a logical');
        if logy
            y = logspace(log10(1),log10(size(z,1)),size(z,1));
        else
            y = 1:size(z,1);
        end        
    end
end
assert(length(x) == size(z,2),'length(x) must equal size(z,2)');
assert(length(y) == size(z,1),'length(y) must equal size(z,1)');
assert(islogical(logx),'logx must be a logical');
assert(islogical(logy),'logy must be a logical');
z = padarray(z,[1 1],'replicate','post');
if logx
    newx = logexpand(x);
else
    newx = linexpand(x);
end
if logy
    newy = logexpand(y);
else
    newy = linexpand(y);
end
p = pcolor(double(real(newx)), double(real(newy)), double(real(z)));
shading flat;
function[ey] = linexpand(y)
%Credit: Benjamin Strom
ey = [y, 0]; 
ey(end) = 2*ey(end-1) - ey(end-2); 
ey = ey - (ey(end-1) - ey(end-2))/2;
function[ex] = logexpand(x)
ex = exp(linexpand(log(x)));
function[x] = prune(x,n)
d = diff(x);
[~,ord] = sort(d);
x(ord(1:n)) = [];