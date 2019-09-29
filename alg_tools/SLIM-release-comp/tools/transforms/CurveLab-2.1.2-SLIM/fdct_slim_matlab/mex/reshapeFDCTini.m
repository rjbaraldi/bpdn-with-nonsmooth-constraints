function [CN, HDR] = reshapeFDCTini(In)

CN = 0;
HDR = cell(size(In));
for i = 1:size(In,2)
  HDR{i} = cell(size(In{i}));
  for j = 1:size(In{i},2)
    tmp = size(In{i}{j});
    HDR{i}{j} = tmp;
    CN = CN + prod(tmp);
  end
end