import pathlib
import pyproj
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import rioxarray as rxr
import streamlit as st
from rasterio import MemoryFile
from owslib.wms import WebMapService
from io import BytesIO
import plotly.graph_objects as go
import plotly.express as px

#matplotlib.use("agg")

CRS_LIST = pyproj.database.query_crs_info()
CRS_STR_LIST = [f"{crs.auth_name}:{crs.code} - {crs.name}" for crs in CRS_LIST]
CRS_DICT = {f"{crs.auth_name}:{crs.code} - {crs.name}": crs for crs in CRS_LIST}
LIDAR_SERVICE_URL=r"https://data.isgs.illinois.edu/arcgis/services/Elevation/IL_Statewide_Lidar_DEM_WGS/ImageServer/WMSServer?request=GetCapabilities&service=WMS"
RASTER_SRC_DICT = {"ISGS Statewide Lidar":LIDAR_SERVICE_URL,
                   "Other Web Service":"get_service_info",
                   "Raster file":"get_file_name"
                   }

def main():
    st.session_state.main_container = st.container()
    with st.sidebar:
        st.title("Raster Sampling on the Web")

        st.selectbox(label="Select raster/elevation data source",
                     options=["ISGS Statewide Lidar", "Other Web Service", "Raster file"],
                     index=0, key='raster_source_select',
                     help="Select the source you would like to use for elevation.",
                     on_change=on_raster_source_change,
                     disabled=True,
                     )

        st.selectbox(label="Raster CRS",
                     options=CRS_STR_LIST,
                     index=CRS_STR_LIST.index("EPSG:3857 - WGS 84 / Pseudo-Mercator"),
                     key='raster_crs',
                     disabled=True)


        raster_expander = st.expander(label='Raster Info.')
        
        st.segmented_control(label='Select point type',
                             options=["Enter coords.",
                                      "CSV File"],
                             selection_mode='single',
                             default="Enter coords.",
                             key="points_source",
                             )

        xcoordCol, ycoordCol = st.columns(spec=[0.5, 0.5], gap='small', border=False)
    
        xcoordCol.text_input(label='X Coordinate',
                      on_change=on_xcoord_change,
                      placeholder="X Coord Value",
                      help="If you enter two comma or space/tab separated values, it will automatically parse the 2nd as your ycoord.",
                      key='xcoord'
                      )

        ycoordCol.text_input(label='Y Coordinate',
                      placeholder="Y Coord Value",
                      key='ycoord'
                      )        
        
        st.selectbox(label="Point CRS",
                     options=CRS_STR_LIST,
                     index=CRS_STR_LIST.index("EPSG:6345 - NAD83(2011) / UTM zone 16N"),
                     key='point_crs')

        point_expander = st.expander(label='Point Info.')

        st.button(label="Sample Elevation Data",
                  key="sample_data_button",
                  icon=":material/background_dot_small:",
                  on_click=get_elevation,
                  type='primary',
                  )

    if hasattr(st.session_state, 'elev_fig'):
        st.session_state.elev_fig

# FUNCTIONS
def on_raster_source_change():
    return

def on_xcoord_change():
    sepList = [',', ' ', '\t']
    if hasattr(st.session_state, 'xcoord'):
        st.session_state.xcoord = st.session_state.xcoord.strip()
        for sep in sepList:
            if sep in st.session_state.xcoord:
                xCoordList = str(st.session_state.xcoord).split(sep)

                if len(xCoordList) > 1:
                    st.session_state.xcoord = xCoordList[0].strip()
                    st.session_state.ycoord = xCoordList[1].strip()
                break

def get_elevation(coords=None, 
                  elevation_col_name='elevation', 
                  xcoord_col_name='xcoord', ycoord_col_name='ycoord', 
                  points_crs=None, 
                  elevation_source=None, elev_source_type='service', 
                  raster_crs=None, show_plot=True):

    if coords is None:
        if st.session_state.points_source=='Enter coords.':
            coords = (st.session_state.xcoord, st.session_state.ycoord)
    

    # Leave everything below here alone
    if elevation_source is None:
        elevation_source = RASTER_SRC_DICT[st.session_state['raster_source_select']]

    if points_crs is None:
        points_crs = int(CRS_DICT[st.session_state.point_crs].code)

    if raster_crs is None:
        raster_crs = int(CRS_DICT[st.session_state.raster_crs].code)

    elev_wms = wms_statewide_lidar_url = elevation_source

    coordTransformer = pyproj.Transformer.from_crs(points_crs, raster_crs)
    wms_statewide_lidar_url = elev_wms

    if isinstance(coords, (tuple, list)):
        xcoord, ycoord = coords
        xcoord_WMS, ycoord_WMS = coordTransformer.transform(xcoord, ycoord)
        minX = maxX = xcoord_WMS
        minY = maxY = ycoord_WMS
    elif isinstance(coords, (pd.DataFrame, gpd.GeoDataFrame, str, pathlib.Path)):
        if isinstance(coords, (pathlib.Path, str)):
            coords = pd.read_csv(coords)
        xcoord = coords[xcoord_col_name]
        ycoord = coords[ycoord_col_name]
        xcoord_WMS, ycoord_WMS = coordTransformer.transform(xcoord, ycoord)

        minX = min(xcoord_WMS)
        minY = min(ycoord_WMS)
        maxX = max(xcoord_WMS)
        maxY = max(ycoord_WMS)


    xPad = (maxX-minX)*0.1
    yPad = (maxY-minY)*0.1

    if xPad == 0:
        xPad = maxX * 0.01
        if abs(xPad) > 1000:
            xPad = 1000
    if yPad == 0:
        yPad = maxY * 0.01
        if abs(yPad) > 1000:
            yPad = 1000
        
    if elev_source_type=='service':
        wms = WebMapService(wms_statewide_lidar_url)

        layer_name = '0'
        layer = wms[layer_name]

        bbox = (minX-xPad, minY-yPad, maxX+xPad, maxY+yPad)
        img = wms.getmap(
            layers=['IL_Statewide_Lidar_DEM_WGS:None'],
            srs='EPSG:3857',
            bbox=bbox,
            size=(256, 256),
            format='image/tiff',
            transparent=True
            )

        try:
            bio = BytesIO(img.read())
            lidarData_rxr = rxr.open_rasterio(bio)
            lidarData = lidarData_rxr.rio.reproject(points_crs)
        except Exception as e:
            print(f"BytesIO method failed: {e}")

        #with MemoryFile(img) as memfile:
            #with memfile.open() as dataset:
        #        lidarData_rxr = rxr.open_rasterio(dataset)
        #lidarData = lidarData_rxr.rio.reproject(points_crs)
    elif elev_source_type == 'file':
        lidarData_rxr = rxr.open_rasterio(elevation_source)

        lidarData = lidarData_rxr.sel(x=slice(minX-xPad, maxX+xPad), y=slice(maxY+yPad, minY-yPad))
        lidarData = lidarData.rio.reproject(points_crs)
    
    minLidarVal = lidarData.min().values
    maxLidarVal = lidarData.max().values
    lidarValRange = maxLidarVal - minLidarVal
    
    vMin = minLidarVal + 0.2*lidarValRange
    vMax = maxLidarVal - 0.2*lidarValRange
    
    if isinstance(coords, (tuple, list)):
        elev_ft = float(lidarData.sel(x=xcoord, y=ycoord, method='nearest').values[0])
        elev_m = elev_ft * 0.3048
        st.session_state.main_container.text(f"elev {elev_ft:.2f} ft    | {elev_m:.2f}m")

        if show_plot:
            #fig, ax = plt.subplots()#figsize=(1,1))
            #lidarData[0].plot(ax=ax, alpha=0.75, cmap='gist_earth', vmin=vMin, vmax=vMax)
            data = lidarData[0].values
            x_coords = lidarData.x.values  
            y_coords = lidarData.y.values

            # Create heatmap
            fig = go.Figure(data=go.Heatmap(
                z=data,
                x=x_coords,
                y=y_coords,
                colorscale='Geyser',
                zmin=vMin,
                zmax=vMax,
                colorbar=dict(title="Elevation (ft)")
            ))
            
            # Add the point marker
            fig.add_trace(go.Scatter(
                x=[xcoord],
                y=[ycoord],
                mode='markers+text',
                marker=dict(
                    symbol='star',
                    size=15,
                    color='red',
                    line=dict(width=2, color='black')
                ),
                text=[f"{elev_ft:.1f} ft<br>{elev_m:.2f} m"],
                textposition="top right",
                textfont=dict(color="red", size=12),
                name=f"Point Elevation",
                showlegend=True
            ))
            
            # Update layout
            fig.update_layout(
                title="Lidar Elevation Data",
                xaxis_title="X Coordinate",
                yaxis_title="Y Coordinate",
                width=700,
                height=500,
                showlegend=True
            )
            
            #st.session_state.main_container.pyplot(fig)
            st.session_state.main_container.plotly_chart(fig)

        coords = (coords[0], coords[1], round(elev_m, 3))
        
    elif isinstance(coords, (pd.DataFrame, gpd.GeoDataFrame)):
        elev_ft = []
        for i, x in enumerate(xcoord):
            y = ycoord[i]
            elev_ft.append(lidarData.sel(x=x, y=y, method='nearest').values[0])
        elev_ft = pd.DataFrame(elev_ft)
        
        elev_m = elev_ft * 0.3048
        coords[elevation_col_name] = elev_m
        
        if show_plot:
            fig, ax = plt.subplots()
            lidarData[0].plot(ax=ax, alpha=0.75, cmap='gist_earth', vmin=vMin, vmax=vMax)
            ax.scatter(x=coords[xcoord_col_name], y=coords[ycoord_col_name], c='k', marker='*',)
            #with st.session_state.main_container:
            #    st.pyplot(fig)

    st.session_state.elev_fig = fig

    return #fig, coords



if __name__ == "__main__":
    main()