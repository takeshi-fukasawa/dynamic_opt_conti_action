function [] = progress(i)
% report progress
  if mod(i, 50) == 0;
    % "Computed: ";;compact(i);;"\r";;
    disp(sprintf('  Computed: %d', i));
    end
end
