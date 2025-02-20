function [out1] = acronym(et,it)
  disp(['Model: ' et ', investment in ' it]);
  if strcmp(it, 'QUALITY')
    if strcmp(et, 'COMPETITION')  out1 = 'cb';
    elseif strcmp(et, 'MONOPOLY') out1 = 'mb';
    elseif strcmp(et, 'PLANNER')  out1 = 'sb';
      end
  elseif strcmp(it, 'COST')
    if strcmp(et, 'COMPETITION')  out1 = 'cc';
    elseif strcmp(et, 'MONOPOLY') out1 = 'mc';
    elseif strcmp(et, 'PLANNER')  out1 = 'sc';
      end
  elseif strcmp(it, 'CAPACITY')
    if strcmp(et, 'COMPETITION')  out1 = 'cp';
    elseif strcmp(et, 'MONOPOLY') out1 = 'mp';
    elseif strcmp(et, 'PLANNER')  out1 = 'sp';
      end
    end

