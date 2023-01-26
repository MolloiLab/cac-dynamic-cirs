### A Pluto.jl notebook ###
# v0.19.18

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 0f1ab7a0-98ec-11ed-0e70-356525ed2545
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using Revise, PlutoUI, DICOM, DICOMUtils, OrthancTools, CairoMakie, LinearAlgebra, Images, StatsBase, Unitful, ActiveContours
end

# ╔═╡ 36c45072-0063-4c2d-ab47-201103f62599
TableOfContents()

# ╔═╡ 8070c932-5d17-4e46-bd2f-566693774599
md"""
# Load DICOMs
"""

# ╔═╡ 7046e8fc-dd3a-433c-8ef2-abb3518b3ce9
ip_address = "128.200.49.26"

# ╔═╡ 15a3fe24-319e-43f8-8ff9-0203ea9b4375
accession_number = "2475"

# ╔═╡ 78f654ed-1f1f-4141-95f3-1b03bb7bf051
series_num = "4"

# ╔═╡ 18eba487-9d8a-4a56-829f-0146692e4ab4
instance_number = 1

# ╔═╡ a2fb46b1-a881-4f14-94d6-2729177f6d37
output_dir = "/Users/daleblack/Documents/dcm_dir"

# ╔═╡ fbc69139-9839-4ef9-bbff-b6d18bc87b72
dcms = dcmdir_parse(output_dir)

# ╔═╡ 7800826f-a9a6-491f-bbe2-cd422a43f6a9
dcm_arr = load_dcm_array(dcms);

# ╔═╡ fe16aad5-65f1-4a0a-ad0b-9ab4eb55e1bb
header = dcms[1].meta

# ╔═╡ fff1e879-3f96-423d-866b-aa3458a8989e
get_pixel_size(header)

# ╔═╡ e0766843-59ea-4a83-8de4-c5ce6952f246
md"""
## Visualize
"""

# ╔═╡ 57b1d4a7-5359-496a-9c83-319476d529ec
@bind a PlutoUI.Slider(axes(dcm_arr, 3), default=130, show_value=true)

# ╔═╡ bdc601d2-7f60-4dbd-af8b-1bdba23e9a61
let
	f = Figure()

	CairoMakie.Axis(f[1, 1])
	heatmap!(dcm_arr[:, :, a], colormap=:grays)

	f
end

# ╔═╡ 6d62bb43-80a3-4e3f-a89a-6282957b17ab
md"""
# Segment Heart
"""

# ╔═╡ 5b76ce9a-a501-49ea-8f66-ac2b53be7370
md"""
## Filter Step

Segment the initial mask used for the Chan Vese level set by filtering the scan for voxels between 30 and 60 HU. Use this to find the center of the heart tissue and then mask a circle about twice the radius of the true heart. This mask then gets used in the Chan Vese level set.
"""

# ╔═╡ be18a03a-09eb-4155-bacb-4fa10556e7ce
function segment_heart(dcm_array, limits=(40, 70))
	# Find the indices of the elements in the array that are between 30 and 60
	selected_inds = findall((limits[1] .<= dcm_array) .& (dcm_array.<= limits[2]))

	# Create boolean array from cartesian indices
	bool_arr = zeros(size(dcm_array))
	for i in selected_inds
		bool_arr[i] = 1
	end

	# Use connected component labeling to identify and label all connected components
	cc_labels = label_components(bool_arr)

	# Use the countmap function to count the number of occurrences of each value in the array, excluding 0
	counts = countmap(cc_labels[cc_labels .!= 0])
	
	# # Find the value with the most occurrences
	most_common_value, _ = sort(collect(pairs(counts)), by=x->x[2], rev=true)
	
	# # Find the indices of the most common value in the original array
	most_common_indices = findall(cc_labels .== most_common_value[1])

	# Create boolean array from new cartesian indices
	bool_arr2 = zeros(size(dcm_array))
	for i in most_common_indices
		bool_arr2[i] = 1
	end
	centroids = Int.(round.(component_centroids(label_components(bool_arr2))[2]))
	
	return bool_arr2, centroids
end

# ╔═╡ 6234f1e4-1799-419b-b345-25803e948316
function create_circle_mask(img, centroids, radius)
    width, height = size(img)
    # initialize level set function with all zeros
    level_set = zeros(size(img))

    # define the center and radius of the initial curve
    center_x, center_y = centroids[1], centroids[2]

    # set all pixels inside the initial curve to -1 and all pixels outside the curve to 1
    for x in 1:width, y in 1:height
        if ((x - center_x)^2 + (y - center_y)^2) <= radius^2
            level_set[x, y] = 1
        else
            level_set[x, y] = 0
        end
    end
    Bool.(level_set)
end

# ╔═╡ 6df5166e-0056-49d6-aaa3-90134ddead90
bool_arr, centroids = segment_heart(dcm_arr, (40, 70));

# ╔═╡ 1b06d20a-7781-4f64-b6c1-1a28ee61463d
dcm_slice = dcm_arr[:, :, 10];

# ╔═╡ 2f479f9f-5e97-4fa2-98b0-ee6ffe225a88
circle_mask = create_circle_mask(dcm_slice, centroids[1:2], 130);

# ╔═╡ 4ff84d02-d4de-4662-aa5c-2049eebccd44
md"""
## Chan Vese

Using the rough mask from above, filter out the checkerboard level set to only segment within this mask. Then erode the segmentation so that the largest connected component is likely the center circle. Use this to find the center of that circle
"""

# ╔═╡ ed458dbc-39d6-48d2-b86e-c128cc8c1b5a
function erode_mask(img, num_erosions)
	new_img = img
	i = 0
	while i < num_erosions
		new_img = erode(new_img)
		i += 1
	end
	return new_img
end

# ╔═╡ 6b2b986d-04fa-4b7c-b239-b1ce277de9a1
begin
	checkerboard = init_checkerboard(size(dcm_slice), 5);
	checkerboard = circle_mask .* checkerboard;
	for i in axes(checkerboard, 1)
		for j in axes(checkerboard, 2)
			if checkerboard[i, j] == 0
				checkerboard[i, j] = -1
			end
		end
	end
end;	

# ╔═╡ c23066e3-a24b-4cf0-9c1c-116d7472ce0b
segmentation, phi, energies = chan_vese(dcm_slice, checkerboard; mu=0.25, max_iter=200);

# ╔═╡ 4ffb82a7-aa38-4e3d-83d6-9ecc3973b506
let
	f = Figure()

	CairoMakie.Axis(
		f[1, 1],
		title="Original Image"
	)
	heatmap!(dcm_slice, colormap=:grays)
	
	CairoMakie.Axis(
		f[1, 2],
		title="Original Level Set"
	)
	heatmap!(checkerboard, colormap=:grays)


	CairoMakie.Axis(
		f[2, 1],
		title="Final Level Set"
	)
	heatmap!(phi, colormap=:grays)
	
	CairoMakie.Axis(
		f[2, 2],
		title="Chan-Vese Segmentation"
	)
	heatmap!(segmentation, colormap=:grays)
	
	f
end

# ╔═╡ 3bb934e6-4e09-4dee-87f7-b7e1df749fb6
begin
	segmentation_eroded = segmentation
	while true
		segmentation_eroded = erode(segmentation_eroded)
		cc = label_components(segmentation_eroded)
		if length(unique(cc)) <= 2
			break
		end
	end
	
	segmentation_eroded_lbl = label_components(segmentation_eroded)
	
	centroids2 = component_centroids(segmentation_eroded_lbl)[2]
end

# ╔═╡ ddd3a997-ebe2-42fc-bb7a-c594fabf770f
component_centroids(segmentation_eroded_lbl)[2]

# ╔═╡ 62c2f9ae-1b58-43ba-9320-318b7b020c9a
begin
	heart_mask = create_circle_mask(dcm_slice, centroids2, 73)
	idxs = findall(isone, heart_mask)
	idxs = getindex.(idxs, [1 2])
end;

# ╔═╡ 70c98c2b-c281-4a00-8dc1-a2ff330bcb74
heatmap(segmentation_eroded)

# ╔═╡ d6ec0afe-afbd-4667-b179-466439864569
md"""
## Visualize
"""

# ╔═╡ 673afc8d-21c1-4592-874c-212c3d704952
@bind b PlutoUI.Slider(axes(dcm_arr, 3), default=130, show_value=true)

# ╔═╡ 00925c19-35e1-4e49-bcd4-2ef0b671a4ee
let
	f = Figure()

	CairoMakie.Axis(f[1, 1])
	heatmap!(dcm_arr[:, :, b], colormap=:grays)
	scatter!(idxs; markersize=2, color=(:red, 0.5))

	f
end

# ╔═╡ Cell order:
# ╠═0f1ab7a0-98ec-11ed-0e70-356525ed2545
# ╠═36c45072-0063-4c2d-ab47-201103f62599
# ╟─8070c932-5d17-4e46-bd2f-566693774599
# ╠═7046e8fc-dd3a-433c-8ef2-abb3518b3ce9
# ╠═15a3fe24-319e-43f8-8ff9-0203ea9b4375
# ╠═78f654ed-1f1f-4141-95f3-1b03bb7bf051
# ╠═18eba487-9d8a-4a56-829f-0146692e4ab4
# ╠═a2fb46b1-a881-4f14-94d6-2729177f6d37
# ╠═fbc69139-9839-4ef9-bbff-b6d18bc87b72
# ╠═7800826f-a9a6-491f-bbe2-cd422a43f6a9
# ╠═fe16aad5-65f1-4a0a-ad0b-9ab4eb55e1bb
# ╠═fff1e879-3f96-423d-866b-aa3458a8989e
# ╟─e0766843-59ea-4a83-8de4-c5ce6952f246
# ╟─57b1d4a7-5359-496a-9c83-319476d529ec
# ╟─bdc601d2-7f60-4dbd-af8b-1bdba23e9a61
# ╟─6d62bb43-80a3-4e3f-a89a-6282957b17ab
# ╟─5b76ce9a-a501-49ea-8f66-ac2b53be7370
# ╠═be18a03a-09eb-4155-bacb-4fa10556e7ce
# ╠═6234f1e4-1799-419b-b345-25803e948316
# ╠═2f479f9f-5e97-4fa2-98b0-ee6ffe225a88
# ╠═6df5166e-0056-49d6-aaa3-90134ddead90
# ╠═1b06d20a-7781-4f64-b6c1-1a28ee61463d
# ╟─4ff84d02-d4de-4662-aa5c-2049eebccd44
# ╟─ed458dbc-39d6-48d2-b86e-c128cc8c1b5a
# ╠═6b2b986d-04fa-4b7c-b239-b1ce277de9a1
# ╠═c23066e3-a24b-4cf0-9c1c-116d7472ce0b
# ╟─4ffb82a7-aa38-4e3d-83d6-9ecc3973b506
# ╠═3bb934e6-4e09-4dee-87f7-b7e1df749fb6
# ╠═ddd3a997-ebe2-42fc-bb7a-c594fabf770f
# ╠═62c2f9ae-1b58-43ba-9320-318b7b020c9a
# ╠═70c98c2b-c281-4a00-8dc1-a2ff330bcb74
# ╟─d6ec0afe-afbd-4667-b179-466439864569
# ╟─673afc8d-21c1-4592-874c-212c3d704952
# ╟─00925c19-35e1-4e49-bcd4-2ef0b671a4ee
