# EOSets

This library aims to simplify any process working with pairs of EO data handled by [EOReader](https://github.com/sertit/eoreader).

Any pair works with a reference product, followed by children products that will align on its grid. 
The pairs will be: [reference, child_1], [reference, child_2], ...

It is assumed that pairs share a common extent where the pair processed will be applied on.
First the extent will be computed, then the reference will be cropped to the wanted outline at the wanted reference, then the children too, and will also be collocated to the reference product.

This library aims to implement several processes, such as:
- pansharpening for multi-product data such as DIMAP (one product is PAN, the other MS)
- constructing pre/post disaster sets
- constructing common grids between products on which to work
- coregistring products
