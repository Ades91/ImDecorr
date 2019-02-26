function val = clamp(val,minval,maxval)

if isempty(minval) % no lower bound
    map = val > maxval;
	val(map) = maxval;
end

if isempty(maxval) % no upper bound
	map = val < minval;
	val(map) = minval;
end

if ~isempty(minval) && ~isempty(maxval)
    map = val > maxval;
    val(map) = maxval;
    map = val < minval;
    val(map) = minval;
end