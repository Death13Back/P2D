function resultingMatrix = blockDiagonalMatrix(param,varargin)
% blockDiagonalMatrix 提供了创建块对角矩阵的接口


  prev_class  = class(varargin{1});
  
  % 检查为构建块对角线提供的矩阵是否都是同一类的变量
  for i=1:nargin-1
    if(isa(varargin{i},prev_class))
      prev_class=class(varargin{i});
    else
      error('All the inputs arguments must be of the same class')
    end
  end
  
  % 使用符号变量运行模型时，如果使用 Octave，则必须以不同的方式创建块对角线
  if(isa(varargin{1},'casadi.SX') || isa(varargin{1},'casadi.MX'))
    if(param.isMatlab)
      resultingMatrix = blkdiag(varargin{:});
    else
      size_tot                                                          = [0,0];  
      for i=1:nargin-1
        size_tot = size_tot + size(varargin{i});
      end
      resultingMatrix                                                             = SX.zeros(size_tot);
      start_position = [1,1];
      % 创建块对角线
      for i=1:nargin-1
        resultingMatrix(start_position(1):start_position(1)+size(varargin{i},1)-1,start_position(2):start_position(2)+size(varargin{i},2)-1) = varargin{i};
        start_position(1) = start_position(1) + size(varargin{i},1);
        start_position(2) = start_position(2) + size(varargin{i},2);
      end
    end
  else
    resultingMatrix = blkdiag(varargin{:});
  end
    
end