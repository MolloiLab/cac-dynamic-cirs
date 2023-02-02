### A Pluto.jl notebook ###
# v0.19.22

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
	using Pkg; Pkg.activate("."); Pkg.instantiate()
	using Revise, PlutoUI, DICOM, OrthancTools, CairoMakie, LinearAlgebra, Images, StatsBase, Unitful, Meshes, Interpolations
	using DICOMUtils, ActiveContours
end

# ╔═╡ 36c45072-0063-4c2d-ab47-201103f62599
TableOfContents()

# ╔═╡ 8070c932-5d17-4e46-bd2f-566693774599
md"""
# Load DICOMs
"""

# ╔═╡ 03323a12-99dc-4d49-bf32-6bfd5f9ef102
md"""
## Locate on Orthanc
"""

# ╔═╡ 7046e8fc-dd3a-433c-8ef2-abb3518b3ce9
# ip_address = "128.200.49.26"

# ╔═╡ 62c3ce5f-2558-48e3-b4e1-fe11d300304f
# studies_dict = get_all_studies(ip_address)

# ╔═╡ 15a3fe24-319e-43f8-8ff9-0203ea9b4375
# accession_number = "2475"

# ╔═╡ 12e74893-e3b2-4f5c-8e3f-e51a7e563a4e
# series_dict = get_all_series(studies_dict, accession_number, ip_address)

# ╔═╡ 78f654ed-1f1f-4141-95f3-1b03bb7bf051
# series_num = "2"

# ╔═╡ 59c1bd7f-62db-4736-bf11-fe071aee81d0
# instances_dict = get_all_instances(series_dict, series_num, ip_address)

# ╔═╡ b1f98e1a-b405-4a72-bba9-17165ac00eaa
md"""
## Download from Orthanc
"""

# ╔═╡ a2fb46b1-a881-4f14-94d6-2729177f6d37
output_dir = "/Users/daleblack/Documents/dcm_dir"

# ╔═╡ 18eba487-9d8a-4a56-829f-0146692e4ab4
# instance_number = 1

# ╔═╡ dc4ea23d-08ce-4607-8bcc-96a4358ba2b6
# download_instances(instances_dict, instance_number, output_dir, ip_address)

# ╔═╡ a7771fa7-d4ea-4c41-9ba2-e03f62fcc2ef
md"""
## Load DICOMs
"""

# ╔═╡ fbc69139-9839-4ef9-bbff-b6d18bc87b72
dcms = dcmdir_parse(output_dir)

# ╔═╡ 7800826f-a9a6-491f-bbe2-cd422a43f6a9
dcm_arr = load_dcm_array(dcms);

# ╔═╡ fe16aad5-65f1-4a0a-ad0b-9ab4eb55e1bb
header = dcms[1].meta;

# ╔═╡ fff1e879-3f96-423d-866b-aa3458a8989e
pixel_size = get_pixel_size(header)

# ╔═╡ e0766843-59ea-4a83-8de4-c5ce6952f246
md"""
## Visualize
"""

# ╔═╡ 57b1d4a7-5359-496a-9c83-319476d529ec
@bind a PlutoUI.Slider(axes(dcm_arr, 3), default=160, show_value=true)

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
    # initialize mask with all zeros
    mask = zeros(size(img))

    # define the center of the circle
    x0, y0 = centroids[1], centroids[2]

    # set all pixels inside the circle to 1 and all pixels outside the circle to 0
    for x in axes(img, 1), y in axes(img, 2)
        if ((x - x0)^2 + (y - y0)^2) <= radius^2
            mask[x, y] = 1
        else
            mask[x, y] = 0
        end
    end
    return Bool.(mask)
end

# ╔═╡ 6df5166e-0056-49d6-aaa3-90134ddead90
bool_arr, centroids = segment_heart(dcm_arr, (40, 70));

# ╔═╡ 1b06d20a-7781-4f64-b6c1-1a28ee61463d
dcm_slice = dcm_arr[:, :, centroids[3]];

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

# ╔═╡ 7fd4884d-92c9-46a6-bbc8-49b3f2f67dd7


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
		if length(unique(cc)) == 2
			break
		end
		if (length(unique(cc)) <= 1)
			@warn "Erosion failed; only background is left"
			break
		end
	end
	
	segmentation_eroded_lbl = label_components(segmentation_eroded)
	
	centroids2 = component_centroids(segmentation_eroded_lbl)[2]
end

# ╔═╡ 62c2f9ae-1b58-43ba-9320-318b7b020c9a
begin
	heart_mask = create_circle_mask(dcm_slice, centroids2, 73)
	idxs = findall(isone, heart_mask)
	idxs = getindex.(idxs, [1 2])
end;

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

# ╔═╡ dec18cd0-720c-4a81-9eaf-8f1e254e414f
md"""
# Mask Heart Function
"""

# ╔═╡ 6b350a05-05f7-421c-8732-cd401eb05ea9
function mask_heart(
	dcm_array, radius1, radius2;
	limits=(40, 70), mu=0.25, lambda1=1, lambda2=1, tol=1e-3, max_iter=200, dt=0.5
)

	## -- Filter Heart Tissue --##
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
	dcm_slice = dcm_arr[:, :, centroids[3]]
	circle_mask = create_circle_mask(dcm_slice, centroids[1:2], radius1)

	## -- Chan Vese --##
	checkerboard = init_checkerboard(size(dcm_slice), 5);
	checkerboard = circle_mask .* checkerboard;
	for i in axes(checkerboard, 1)
		for j in axes(checkerboard, 2)
			if checkerboard[i, j] == 0
				checkerboard[i, j] = -1
			end
		end
	end

	segmentation, phi, energies = chan_vese(
		dcm_slice, checkerboard; 
		mu=mu, lambda1=lambda1, lambda2=lambda2, tol=tol, max_iter=max_iter, dt=dt);

	while true
		segmentation = erode(segmentation)
		cc = label_components(segmentation)
		if length(unique(cc)) == 2
			break
		end
		if (length(unique(cc)) <= 1)
			@warn "Erosion failed; only background is left"
			break
		end
	end
	
	segmentation_lbl = label_components(segmentation)
	centroids2 = component_centroids(segmentation_lbl)[2]

	heart_mask = create_circle_mask(dcm_slice, centroids2, radius2)
	return heart_mask
end

# ╔═╡ ef32e76b-39dc-4f23-9db8-d859cbaa7cea
begin
	heart_mask2 = mask_heart(dcm_arr, 130, 70)
	idxs2 = findall(isone, heart_mask2)
	idxs2 = getindex.(idxs2, [1 2])
end;

# ╔═╡ 4cd425c1-e6ef-4bd8-b6fd-c79fdfe7b7fe
md"""
## Visualize
"""

# ╔═╡ 928dc26c-c133-42a5-a18b-1b592fac70b4
@bind c PlutoUI.Slider(axes(dcm_arr, 3), default=130, show_value=true)

# ╔═╡ b6abd768-1c50-4bef-8cac-ae8c5b23939f
let
	f = Figure()

	CairoMakie.Axis(f[1, 1])
	heatmap!(dcm_arr[:, :, c], colormap=:grays)
	scatter!(idxs2; markersize=2, color=(:red, 0.5))

	f
end

# ╔═╡ 0ae5c8d2-cc9a-4795-a0dc-21227bee8547
md"""
# Inserts
"""

# ╔═╡ 3f8d9d36-ca54-41ea-94d9-2a3656dc9622
dcm_heart = dcm_arr .* heart_mask2;

# ╔═╡ 6355709c-c97a-4d13-9337-e413eb309f08
@bind d PlutoUI.Slider(axes(dcm_heart, 3), default=130, show_value=true)

# ╔═╡ cf18c58f-9270-42f8-bd32-a109459171a3
let
	f = Figure()

	CairoMakie.Axis(f[1, 1])
	heatmap!(dcm_heart[:, :, d], colormap=:grays)
	
	f
end

# ╔═╡ 4764a57a-2210-4398-b017-451922913c5c
md"""
## Endpoints
"""

# ╔═╡ 4213eae9-e95d-46e0-a617-e882f893853f
function find_heart_endpoints(dcm_heart, air_threshold=-2000)
	# Find the indices of the elements in the array that are between 30 and 60
	selected_inds = findall(dcm_heart .<= air_threshold)
	
	# Create boolean array from cartesian indices
	local bool_arr4 = zeros(size(dcm_heart))
	for i in selected_inds
		bool_arr4[i] = 1
	end

	local cc_labels 
	while true
		# Use connected component labeling to identify and label all connected components
		cc_labels = label_components(bool_arr4)
		if length(unique(cc_labels)) == 3
			break
		end
		
		bool_arr4 = erode(bool_arr4)
		if (length(unique(cc_labels)) <= 2)
			@warn "Erosion failed; only background is left"
			break
		end
	end
	
	begining_slice = findlast(isone, cc_labels)[3] + 1
	end_slice = findfirst(x->x==2, cc_labels)[3] - 1
	return begining_slice, end_slice
end

# ╔═╡ 14eadcb4-f842-4e96-a399-5246c30c4f44
begining_slice, end_slice = find_heart_endpoints(dcm_heart)

# ╔═╡ 863430eb-ac14-4fa0-96db-b57b0c1df70d
md"""
## Plane Fitting
"""

# ╔═╡ 54069a4c-1a24-4a91-91a1-654846141040
function find_heart_plane(dcm_heart, endpoints, air_threshold=-300)
	# Find the indices of the elements in the array that are less than
	# air_threshold and filter to exclude beginning and end slices
	remove = [collect(1:endpoints[1])..., collect(endpoints[2]:320)...]
	selected_inds = findall(dcm_heart .<= air_threshold)
	for r in remove
		selected_inds = [i for i in selected_inds if i.I[3] != r]
	end
	return selected_inds
end

# ╔═╡ a2c984ce-0922-4ba7-8a7b-af40e1c0cc2b
begin
	idx_plane = find_heart_plane(dcm_heart, (begining_slice, end_slice))
	idx_plane = getindex.(idx_plane, [1 2 3])
end

# ╔═╡ dfe9949d-f336-483f-b437-8667992537d2
@bind e PlutoUI.Slider(axes(dcm_heart, 3), default=125, show_value=true)

# ╔═╡ 02482499-8cc1-45ff-835a-1c4e5f7ee0aa
let
	idx_plane_2d = findall(idx_plane[:, 3] .== e)
	f, l = first(idx_plane_2d), last(idx_plane_2d)
	idx_plane_2d = idx_plane[f:l, 1:2]
	
	f = Figure()
	
	ax = CairoMakie.Axis(f[1, 1])
	heatmap!(dcm_heart[:, :, e], colormap=:grays)
	scatter!(idx_plane_2d; markersize=2, color=:red)

	f
end

# ╔═╡ 4f0aecf0-6a2d-4a18-a191-954109cad93d
# barycenter of the points
# compute centered coordinates
G = sum(idx_plane; dims=1) ./ size(idx_plane, 1)

# ╔═╡ ad04b838-2aad-4f5a-9b48-f184c9b08675
# run SVD
u, s, vh = svd(idx_plane .- G);

# ╔═╡ 77a52824-c99f-4ac8-afb2-53ad131d66a2
unit_normal = vh[3, :] # extract normal vector

# ╔═╡ 6b78e611-279b-40ca-a020-5b165741f799
pt = idx_plane[20, :]

# ╔═╡ 060aa588-290f-44e2-8743-f5c180c265fb
plane = Plane(Meshes.Point(pt), Meshes.Vec(unit_normal))

# ╔═╡ 1de5f162-9d6b-4625-8540-990b145db94e
A_itp = interpolate(dcm_arr, BSpline(Cubic(Interpolations.Line(OnGrid()))));

# ╔═╡ e8013c8e-7c3b-4a60-b00f-94089cfd32d8
extrapolation_bc = 0.0

# ╔═╡ 759c81be-cf22-4d0f-834b-66f509c3aa8f
begin
	r = 512/2
	r2 = 320/2
	locations = [plane(i, j) + Meshes.Vec(unit_normal)*k for i in -r:r, j in -r:r, k in -r2:r2]
end

# ╔═╡ 82dcd0fd-4e6b-4ae8-9d29-770be4ea5ad5
etp = extrapolate(A_itp, Periodic());

# ╔═╡ fa59639a-74b1-4b15-bee3-28f737451560
test_arr = [etp(loc.coords...) for loc in locations]

# ╔═╡ 2dd189ed-6bea-4030-acb6-0d7bbf0f0640
@bind f PlutoUI.Slider(axes(test_arr, 3), default=125, show_value=true)

# ╔═╡ ad9367c0-ecd2-4307-baee-4f2615eb1709
heatmap(test_arr[:, :, f], colormap=:grays)

# ╔═╡ Cell order:
# ╠═0f1ab7a0-98ec-11ed-0e70-356525ed2545
# ╠═36c45072-0063-4c2d-ab47-201103f62599
# ╟─8070c932-5d17-4e46-bd2f-566693774599
# ╟─03323a12-99dc-4d49-bf32-6bfd5f9ef102
# ╠═7046e8fc-dd3a-433c-8ef2-abb3518b3ce9
# ╠═62c3ce5f-2558-48e3-b4e1-fe11d300304f
# ╠═15a3fe24-319e-43f8-8ff9-0203ea9b4375
# ╠═12e74893-e3b2-4f5c-8e3f-e51a7e563a4e
# ╠═78f654ed-1f1f-4141-95f3-1b03bb7bf051
# ╠═59c1bd7f-62db-4736-bf11-fe071aee81d0
# ╟─b1f98e1a-b405-4a72-bba9-17165ac00eaa
# ╠═a2fb46b1-a881-4f14-94d6-2729177f6d37
# ╠═18eba487-9d8a-4a56-829f-0146692e4ab4
# ╠═dc4ea23d-08ce-4607-8bcc-96a4358ba2b6
# ╟─a7771fa7-d4ea-4c41-9ba2-e03f62fcc2ef
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
# ╠═7fd4884d-92c9-46a6-bbc8-49b3f2f67dd7
# ╠═6b2b986d-04fa-4b7c-b239-b1ce277de9a1
# ╠═c23066e3-a24b-4cf0-9c1c-116d7472ce0b
# ╟─4ffb82a7-aa38-4e3d-83d6-9ecc3973b506
# ╠═3bb934e6-4e09-4dee-87f7-b7e1df749fb6
# ╠═62c2f9ae-1b58-43ba-9320-318b7b020c9a
# ╟─d6ec0afe-afbd-4667-b179-466439864569
# ╟─673afc8d-21c1-4592-874c-212c3d704952
# ╟─00925c19-35e1-4e49-bcd4-2ef0b671a4ee
# ╟─dec18cd0-720c-4a81-9eaf-8f1e254e414f
# ╠═6b350a05-05f7-421c-8732-cd401eb05ea9
# ╠═ef32e76b-39dc-4f23-9db8-d859cbaa7cea
# ╟─4cd425c1-e6ef-4bd8-b6fd-c79fdfe7b7fe
# ╟─928dc26c-c133-42a5-a18b-1b592fac70b4
# ╟─b6abd768-1c50-4bef-8cac-ae8c5b23939f
# ╟─0ae5c8d2-cc9a-4795-a0dc-21227bee8547
# ╠═3f8d9d36-ca54-41ea-94d9-2a3656dc9622
# ╟─6355709c-c97a-4d13-9337-e413eb309f08
# ╟─cf18c58f-9270-42f8-bd32-a109459171a3
# ╟─4764a57a-2210-4398-b017-451922913c5c
# ╠═4213eae9-e95d-46e0-a617-e882f893853f
# ╠═14eadcb4-f842-4e96-a399-5246c30c4f44
# ╟─863430eb-ac14-4fa0-96db-b57b0c1df70d
# ╠═54069a4c-1a24-4a91-91a1-654846141040
# ╠═a2c984ce-0922-4ba7-8a7b-af40e1c0cc2b
# ╟─dfe9949d-f336-483f-b437-8667992537d2
# ╟─02482499-8cc1-45ff-835a-1c4e5f7ee0aa
# ╠═4f0aecf0-6a2d-4a18-a191-954109cad93d
# ╠═ad04b838-2aad-4f5a-9b48-f184c9b08675
# ╠═77a52824-c99f-4ac8-afb2-53ad131d66a2
# ╠═6b78e611-279b-40ca-a020-5b165741f799
# ╠═060aa588-290f-44e2-8743-f5c180c265fb
# ╠═1de5f162-9d6b-4625-8540-990b145db94e
# ╠═e8013c8e-7c3b-4a60-b00f-94089cfd32d8
# ╠═759c81be-cf22-4d0f-834b-66f509c3aa8f
# ╠═82dcd0fd-4e6b-4ae8-9d29-770be4ea5ad5
# ╠═fa59639a-74b1-4b15-bee3-28f737451560
# ╟─2dd189ed-6bea-4030-acb6-0d7bbf0f0640
# ╠═ad9367c0-ecd2-4307-baee-4f2615eb1709
